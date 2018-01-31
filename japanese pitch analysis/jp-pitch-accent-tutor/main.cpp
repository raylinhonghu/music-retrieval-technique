#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>

#include <cstdio>
#include <fstream>
#include <codecvt>
#include <atomic>
#include <algorithm>
#include <complex>
#include <cassert>

// Reads a file encoded in UTF-8 with automatic conversion to UTF-16.
// UTF-16 is easier for Japanese than UTF-8. Not perfect, but good enough.
// Note "wide string" == UTF-16 here
std::wifstream OpenInputUTF8FileAsWide(const char* filename)
{
    std::wifstream wif(filename);
    // note: no need to call delete, since the std locale object does it.
    wif.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));
    return std::move(wif);
}

class KatakanaOps
{
    static const wchar_t kFirstKatakana = 0x30A0;
    static const wchar_t kLastKatakana = 0x30FF;

    static const std::string* GetRomajiStringTable()
    {
        static const std::string skTable[] = {
            "=",  "a",   "a",  "i",  "i",   "u",  "u",  "e",   "e",  "o",  "o",  "ka", "ga", "ki", "gi", "ku",
            "gu", "ke",  "ge", "ko", "go",  "sa", "za", "shi", "ji", "su", "zu", "se", "ze", "so", "zo", "ta",
            "da", "chi", "ji", "ltsu",  "tsu", "zu", "te", "de",  "to", "do", "na", "ni", "nu", "ne", "no", "ha",
            "ba", "pa",  "hi", "bi", "pi",  "fu", "bu", "pu",  "he", "be", "pe", "ho", "bo", "po", "ma", "mi",
            "mu", "me",  "mo", "ya", "ya",  "yu", "yu", "yo",  "yo", "ra", "ri", "ru", "re", "ro", "wa", "wa",
            "_",  "_",   "wo", "n",  "v",   "_",  "_",  "_",   "_",  "_",  "vo", " ",  "chouonpu",  "_",  "_",  "_"
        };

        return skTable;
    }

public:
    static int NextMora(const std::wstring& katakana, int currMora)
    {
        assert(currMora < (int)katakana.size());

        static const bool skIsMoraStartTable[] = {
            0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
            1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0
        };

        int next;
        for (next = currMora + 1; next < (int)katakana.size(); next++)
        {
            wchar_t kana = katakana[next];
            assert(kana >= kFirstKatakana && kana <= kLastKatakana);
            if (skIsMoraStartTable[kana - kFirstKatakana])
            {
                break;
            }
        }

        return next;
    }

    static int MoraCount(const std::wstring& katakana)
    {
        int cnt = 0;
        for (int currMora = 0; currMora < (int)katakana.size(); currMora = NextMora(katakana, currMora))
        {
            cnt++;
        }
        return cnt;
    }

    static std::string KatakanaToRomaji(const std::wstring& katakana)
    {
        const std::string* romajiTable = GetRomajiStringTable();

        std::string s;
        for (int currMora = 0; currMora < (int)katakana.size();)
        {
            int nextMora = NextMora(katakana, currMora);

            assert(nextMora - currMora == 1 || nextMora - currMora == 2);

            if (nextMora - currMora == 1)
            {
                s += romajiTable[katakana[currMora] - kFirstKatakana];
            }
            else
            {
                s += romajiTable[katakana[currMora] - kFirstKatakana].substr(0, 1)
                    + romajiTable[katakana[currMora + 1] - kFirstKatakana];
            }

            currMora = nextMora;
        }

        return s;
    }

};

struct PitchDictionaryEntry
{
    std::wstring kanji;
    std::wstring katakana;
    std::string romaji;
    std::vector<uint8_t> pitchcode;
};

class PitchAccentDictionary
{
    std::vector<PitchDictionaryEntry> m_dictionary;

public:
    bool loadFromFile(const char* filename)
    {
        m_dictionary.clear();

        std::wifstream dictionary = OpenInputUTF8FileAsWide(filename);
        if (!dictionary)
        {
            printf("Failed to open dictionary file %s\n", filename);
            return false;
        }

        PitchDictionaryEntry currEntry;

        std::wstring line;
        while (std::getline(dictionary, line))
        {
            size_t currTokenStart = 0, currTokenEnd = 0;

            // extract token for the kanji of the word
            currTokenEnd = line.find('\t', currTokenStart);
            if (currTokenEnd == std::wstring::npos)
            {
                continue;
            }
            currEntry.kanji.assign(line.c_str() + currTokenStart, line.c_str() + currTokenEnd);
            currTokenStart = currTokenEnd + 1;

            // extract token for the katakana of the word
            currTokenEnd = line.find('\t', currTokenStart);
            if (currTokenEnd == std::wstring::npos)
            {
                continue;
            }
            currEntry.katakana.assign(line.c_str() + currTokenStart, line.c_str() + currTokenEnd);
            currTokenStart = currTokenEnd + 1;

            currEntry.romaji = KatakanaOps::KatakanaToRomaji(currEntry.katakana);

            // extract the pitch code
            currEntry.pitchcode.clear();
            for (size_t i = currTokenStart; i < line.size(); i++)
            {
                wchar_t charcode = line[i];
                
                int intcode = -1;
                if (charcode >= L'0' && charcode <= L'9')
                {
                    intcode = charcode - L'0';
                }
                else if (charcode >= L'A' && charcode <= L'F')
                {
                    intcode = 0xA + charcode - L'A';
                }
                else
                {
                    printf("Invalid pitch code\n");
                    return false;
                }

                currEntry.pitchcode.push_back(intcode);
            }

            m_dictionary.push_back(currEntry);
        }

        return true;
    }

    const PitchDictionaryEntry& at(size_t i) const
    {
        return m_dictionary.at(i);
    }

    size_t size() const
    {
        return m_dictionary.size();
    }
};

// Adapted from https://stackoverflow.com/a/37729648/2752296
class FFT
{
    static constexpr float kPI = 3.14159265359f;

    // pre-allocated up front
    std::vector<std::complex<float>> m_W;
    std::vector<std::complex<float>> m_f2;

    static int log2(int N)    /*function to calculate the log2(.) of int numbers*/
    {
        int k = N, i = 0;
        while (k) {
            k >>= 1;
            i++;
        }
        return i - 1;
    }

    static int reverse(int N, int n)    //calculating revers number
    {
        int j, p = 0;
        for (j = 1; j <= log2(N); j++) {
            if (n & (1 << (log2(N) - j)))
                p |= 1 << (j - 1);
        }
        return p;
    }

    void ordina(std::complex<float>* f1, int N) //using the reverse order in the array
    {
        for (int i = 0; i < N; i++)
            m_f2[i] = f1[reverse(N, i)];
        for (int j = 0; j < N; j++)
            f1[j] = m_f2[j];
    }

    void transform(std::complex<float>* f, int N) //
    {
        ordina(f, N);    //first: reverse order
        m_W[1] = std::polar(1., -2. * kPI / N);
        m_W[0] = 1;
        for (int i = 2; i < N / 2; i++)
            m_W[i] = pow(m_W[1], i);
        int n = 1;
        int a = N / 2;
        for (int j = 0; j < log2(N); j++) {
            for (int i = 0; i < N; i++) {
                if (!(i & n)) {
                    std::complex<float> temp = f[i];
                    std::complex<float> Temp = m_W[(i * a) % (n * a)] * f[i + n];
                    f[i] = temp + Temp;
                    f[i + n] = temp - Temp;
                }
            }
            n *= 2;
            a = a / 2;
        }
    }

public:
    FFT(int maxN)
        : m_W(maxN)
        , m_f2(maxN)
    { }

    void operator()(std::complex<float>* f, int N, float d)
    {
        transform(f, N);
        for (int i = 0; i < N; i++)
            f[i] *= d; //multiplying by step
    }
};

// According to https://en.wikipedia.org/wiki/Pitch_detection_algorithm#Fundamental_frequency_of_speech
const float kSpeechLowestFrequency = 40.0f;
const float kSpeechHighestFrequency = 600.0f;

class VoiceRecorder : public sf::SoundRecorder
{
    static const int kMaxWindowSize = 2048;
    static_assert(((kMaxWindowSize - 1) & kMaxWindowSize) == 0, "Window size must be a power of two");

    FFT m_fft{ kMaxWindowSize };

    std::vector<std::complex<float>> m_fftWindow;
    std::vector<float> m_windowF0s;

    std::atomic<float> m_currVolume = 0.0f;
    std::atomic<float> m_currPitch = 0.0f;

public:
    VoiceRecorder()
    {
        // allocate window buffer up front
        m_fftWindow.resize(kMaxWindowSize);
    }

    bool onProcessSamples(const sf::Int16* samples, size_t sampleCount) override
    {
        int N = int(sampleCount);

        // calculate average volume
        {
            int64_t avg = 0;
            for (size_t i = 0; i < sampleCount; i++)
            {
                avg += abs(int32_t(samples[i]));
            }
            avg /= sampleCount;

            float favg = avg / float(-int32_t(SHRT_MIN));
            favg = favg < 0.0f ? 0.0f : favg > 1.0f ? 1.0f : favg;

            m_currVolume = favg;
        }

        int numWindows = (N + kMaxWindowSize - 1) / kMaxWindowSize;

        m_windowF0s.resize(numWindows);

        for (int windowIdx = 0; windowIdx < numWindows; windowIdx++)
        {
            int windowSize = std::min(kMaxWindowSize, N - windowIdx * kMaxWindowSize);
            const sf::Int16* window = samples + windowIdx * kMaxWindowSize;

            // ignore silence
            bool silent = true;
            for (int i = 0; i < windowSize; i++)
            {
                if (abs(window[i]) > 100)
                {
                    silent = false;
                    break;
                }
            }

            if (silent)
            {
                m_windowF0s[windowIdx] = 0.0f;
                continue;
            }
            
            // convert window to floats
            for (int i = 0; i < windowSize; i++)
            {
                m_fftWindow[i] = float(window[i]) / float(-SHRT_MIN);
            }

            // get amplitude spectrum of fft
            m_fft(m_fftWindow.data(), windowSize, 1.0);

            // convert to hertz
            for (int i = 0; i < windowSize; i++)
            {
                m_fftWindow[i] = abs(m_fftWindow[i].real() / float(windowSize));
            }

            // find the highest peak in the amplitude spectrum
            int peakIdx = 0;
            float peak = abs(m_fftWindow[peakIdx]);
            for (int i = 1; i < windowSize; i++)
            {
                float amp = abs(m_fftWindow[i]);
                if (amp > peak)
                {
                    peakIdx = i;
                    peak = amp;
                }
            }

            if (peak < 0.01f)
            {
                // too quiet
                m_windowF0s[windowIdx] = 0.0f;
                continue;
            }

            // get the frequency that corresponds to the peak
            float f0 = float(getSampleRate()) * float(peakIdx) / 2.0f / float(windowSize);

            // clamp to reasonable range
            if (f0 < kSpeechLowestFrequency)
            {
                f0 = kSpeechLowestFrequency;
            }
            else if (f0 > kSpeechHighestFrequency)
            {
                f0 = kSpeechHighestFrequency;
            }

            m_windowF0s[windowIdx] = f0;
        }

        // average them out
        float f0 = 0.0f;
        for (float wf0 : m_windowF0s)
        {
            f0 += wf0;
        }
        f0 /= float(m_windowF0s.size());
        
        // clamp
        if (f0 < kSpeechLowestFrequency)
        {
            f0 = kSpeechLowestFrequency;
        }
        else if (f0 > kSpeechHighestFrequency)
        {
            f0 = kSpeechHighestFrequency;
        }

        m_currPitch = f0;

        return true;
    }

    float getCurrentPitch() const
    {
        return m_currPitch;
    }

    float getCurrentVolume() const
    {
        return m_currVolume;
    }
};

class PitchGuideRenderer : public sf::Drawable
{
    static constexpr float kPitchCircleOutlineThickness = 10.0f;
    static constexpr float kPitchCircleRadius = 20.0f;
    static constexpr float kPitchCircleLeftBoundary = 250.0f;
    static constexpr float kPitchCircleTopBoundary = 250.0f;
    static constexpr float kPitchCircleHorizontalSpacing = 200.0f;
    static constexpr float kPitchCircleHighToLowHeightDiff = 150.0f;
    
    static constexpr float kSecondsPerMora = 0.2f;

    struct PitchHistorySample
    {
        float time;
        float pitch;
        float volume;
    };

    const sf::Font& m_font;

    std::vector<sf::Text> m_romajiTexts;
    std::vector<sf::Text> m_katakanaTexts;
    sf::Text m_kanjiText;

    sf::Text m_instructionsText;

    sf::Text m_defaultScoreText;
    sf::Text m_scoreText;
    
    std::vector<sf::CircleShape> m_guideCircles;
    sf::VertexArray m_guideCirclePathVerts;

    bool m_hasCurrEntry = false;
    PitchDictionaryEntry m_currEntry;

    sf::CircleShape m_pitchCircle;
    sf::VertexArray m_pitchCirclePathVerts;

    float m_currPitch = 0.0f;
    float m_currVolume = 0.0f;

    static constexpr float kInitialWordTime = -0.2f;

    float m_currWordTime = kInitialWordTime;
    std::vector<PitchHistorySample> m_currWordPitchHistory;
    bool m_playGuideAnimation = false;

    float m_currScore = NAN;

public:
    explicit PitchGuideRenderer(const sf::Font& font)
        : m_font(font)
    {
        m_kanjiText.setFont(m_font);
        m_kanjiText.setCharacterSize(50);

        m_instructionsText.setFont(m_font);
        m_instructionsText.setCharacterSize(36);
        m_instructionsText.setString(
            "Instructions:\n"
            "  Space: Play/Stop\n"
            "  Left/Right: Switch word\n"
        );
        
        m_defaultScoreText.setFont(m_font);
        m_defaultScoreText.setCharacterSize(36);
        m_defaultScoreText.setString("Accuracy: N/A");

        m_scoreText.setFont(m_font);
        m_scoreText.setCharacterSize(36);

        m_guideCirclePathVerts.setPrimitiveType(sf::Triangles);

        m_pitchCirclePathVerts.setPrimitiveType(sf::Triangles);
    }

    void setCurrentWord(const PitchDictionaryEntry& entry)
    {
        int moraCount = KatakanaOps::MoraCount(entry.katakana);
        m_romajiTexts.resize(moraCount);
        m_katakanaTexts.resize(moraCount);
        
        for (int currMora = 0, moraIdx = 0; currMora < (int)entry.katakana.size(); moraIdx++)
        {
            int nextMora = KatakanaOps::NextMora(entry.katakana, currMora);

            m_romajiTexts[moraIdx].setFont(m_font);
            m_romajiTexts[moraIdx].setCharacterSize(50);
            m_romajiTexts[moraIdx].setString(KatakanaOps::KatakanaToRomaji(entry.katakana.substr(currMora, nextMora - currMora)));

            m_katakanaTexts[moraIdx].setFont(m_font);
            m_katakanaTexts[moraIdx].setCharacterSize(50);
            m_katakanaTexts[moraIdx].setString(entry.katakana.substr(currMora, nextMora - currMora));

            currMora = nextMora;
        }

        for (int currMora = 0, moraIdx = 0; currMora < (int)entry.katakana.size(); moraIdx++)
        {
            int nextMora = KatakanaOps::NextMora(entry.katakana, currMora);

            if (m_romajiTexts[moraIdx].getString() == "ltsu")
            {
                m_romajiTexts[moraIdx].setString(m_romajiTexts.at(moraIdx + 1).getString().substring(0, 1));
            }
            if (m_romajiTexts[moraIdx].getString() == "chouonpu")
            {
                sf::String prev = m_romajiTexts.at(moraIdx - 1).getString();
                m_romajiTexts[moraIdx].setString(prev.substring(prev.getSize() - 1, 1));
            }

            currMora = nextMora;
        }

        m_kanjiText.setString(entry.kanji);

        m_currEntry = entry;
        m_hasCurrEntry = true;
        m_currWordTime = kInitialWordTime;
        m_playGuideAnimation = false;
        m_currWordPitchHistory.clear();
    }

    void setCurrentPitch(float pitch)
    {
        m_currPitch = pitch;
    }

    void setCurrentVolume(float volume)
    {
        m_currVolume = volume;
    }

    void update(float deltaTime, sf::RenderTarget& target)
    {
        const sf::View& view = target.getView();

        const float viewLeftRightPadding = view.getSize().x * 0.05f;
        const float viewTopBottomPadding = view.getSize().y * 0.02f;

        m_instructionsText.setPosition(
            view.getSize().x - m_instructionsText.getGlobalBounds().width - (m_instructionsText.getGlobalBounds().left - m_instructionsText.getPosition().x) - viewLeftRightPadding,
            view.getSize().y - m_instructionsText.getGlobalBounds().height - (m_instructionsText.getGlobalBounds().top - m_instructionsText.getPosition().y) - viewTopBottomPadding);

        if (!m_hasCurrEntry)
        {
            return;
        }

        float currWordTimeDuration = float(m_currEntry.pitchcode.size() + 1) * kSecondsPerMora;

        if (m_playGuideAnimation)
        {
            PitchHistorySample newSample{ m_currWordTime, m_currPitch, m_currVolume };

            size_t oldHistorySize = m_currWordPitchHistory.size();
            if (m_currWordPitchHistory.empty() || m_currWordTime - m_currWordPitchHistory.back().time > 1.0f / 60.0f)
            {
                if (m_currWordPitchHistory.empty())
                {
                    m_currWordPitchHistory.push_back(newSample);
                }
                else
                {
                    PitchHistorySample lastSample = m_currWordPitchHistory.back();
                    bool shouldAdd = newSample.pitch != lastSample.pitch || (newSample.time - lastSample.time > kSecondsPerMora);
                    if (shouldAdd)
                    {
                        m_currWordPitchHistory.push_back(newSample);
                    }
                }
            }
            bool addedToPitchHistory = m_currWordPitchHistory.size() != oldHistorySize;

            m_currWordTime += deltaTime;
            if (m_currWordTime > currWordTimeDuration)
            {
                m_currWordTime = currWordTimeDuration;
                m_playGuideAnimation = false;

                if (!addedToPitchHistory)
                {
                    // always duplicate the last one
                    m_currWordPitchHistory.push_back(newSample);
                }

                // calculate score
                enum
                {
                    is_peak = 1,
                    is_trough = 2,
                    is_increase = 4,
                    is_decrease = 8,
                    is_plateau = 16,
                    is_valley = 32,
                };

                float avgPitch = 0.0f;
                float minPitch = INFINITY;
                float maxPitch = -INFINITY;
                float avgVolume = 0.0f;
                float minVolume = INFINITY;
                float maxVolume = -INFINITY;
                for (const auto& hist : m_currWordPitchHistory)
                {
                    avgPitch += hist.pitch;
                    minPitch = std::min(minPitch, hist.pitch);
                    maxPitch = std::max(maxPitch, hist.pitch);
                    avgVolume += hist.volume;
                    minVolume = std::min(minVolume, hist.volume);
                    maxVolume = std::max(maxVolume, hist.volume);
                }
                avgPitch /= float(m_currWordPitchHistory.size());
                avgVolume /= float(m_currWordPitchHistory.size());
                float pitchDiff = maxPitch - minPitch;
                float volumeDiff = maxVolume - minVolume;

                struct Feature
                {
                    int type;
                    float time;
                };
                std::vector<Feature> features;
                for (int i = 1; i < (int)m_currWordPitchHistory.size() - 1; i++)
                {
                    if (m_currWordPitchHistory[i].volume < minVolume + volumeDiff * 0.1f)
                    {
                        // too quiet
                        continue;
                    }

                    float movement = std::fabsf(m_currWordPitchHistory[i].pitch - m_currWordPitchHistory[i - 1].pitch) + std::fabsf(m_currWordPitchHistory[i + 1].pitch - m_currWordPitchHistory[i].pitch);
                    bool flatenuf = movement == 0.0f || movement < pitchDiff / 3.0f;

                    int feature;

                    if (flatenuf)
                    {
                        if (m_currWordPitchHistory[i].pitch > avgPitch)
                        {
                            feature = is_plateau;
                        }
                        else
                        {
                            feature = is_valley;
                        }
                    }
                    else
                    {
                        if (m_currWordPitchHistory[i - 1].pitch <= m_currWordPitchHistory[i].pitch && m_currWordPitchHistory[i].pitch >= m_currWordPitchHistory[i + 1].pitch)
                        {
                            feature = is_peak;
                        }
                        else if (m_currWordPitchHistory[i - 1].pitch >= m_currWordPitchHistory[i].pitch && m_currWordPitchHistory[i].pitch <= m_currWordPitchHistory[i + 1].pitch)
                        {
                            feature = is_trough;
                        }
                        else if (m_currWordPitchHistory[i - 1].pitch <= m_currWordPitchHistory[i].pitch && m_currWordPitchHistory[i].pitch <= m_currWordPitchHistory[i + 1].pitch)
                        {
                            feature = is_increase;
                        }
                        else if (m_currWordPitchHistory[i - 1].pitch >= m_currWordPitchHistory[i].pitch && m_currWordPitchHistory[i].pitch >= m_currWordPitchHistory[i + 1].pitch)
                        {
                            feature = is_decrease;
                        }
                        else
                        {
                            assert(false);
                        }
                    }
                    features.push_back(Feature{ feature, m_currWordPitchHistory[i].time });
                }

                int goodMoras = 0;
                int currFeature = 0;
                for (int i = 0; i < (int)m_currEntry.pitchcode.size(); i++)
                {
                    bool last = i == 0 ? (m_currEntry.pitchcode[i] & 0x1) : (m_currEntry.pitchcode[i - 1] & 0x1);
                    bool curr = (m_currEntry.pitchcode[i] & 0x1);
                    bool next = i == (int)m_currEntry.pitchcode.size() - 1 ? ((m_currEntry.pitchcode[i] & 0x8) ? false : (m_currEntry.pitchcode[i] & 0x1)) : (m_currEntry.pitchcode[i + 1] & 0x1);
                    int pitchbits = (last ? 4 : 0) | (curr ? 2 : 0) | (next ? 1 : 0);

                    static const int matchMatrix[] = {
                        /* lo, lo, lo */ is_trough | is_decrease | is_increase | is_valley, 
                        /* lo, lo, hi */ is_trough | is_increase | is_valley,
                        /* lo, hi, lo */ is_peak   | is_increase, 
                        /* lo, hi, hi */ is_peak   | is_increase, 
                        /* hi, lo, lo */ is_trough | is_decrease, 
                        /* hi, lo, hi */ is_trough | is_decrease, 
                        /* hi, hi, lo */ is_peak   | is_decrease | is_plateau, 
                        /* hi, hi, hi */ is_trough | is_decrease | is_increase | is_plateau,
                    };

                    int matchbits = matchMatrix[pitchbits];

                    while (currFeature < (int)features.size())
                    {
                        bool closeenuf = std::fabsf(features[currFeature].time - (i * kSecondsPerMora - kInitialWordTime)) < kSecondsPerMora / 2.0f;
                        bool currmatched = (features[currFeature].type & matchbits) != 0 && closeenuf;

                        currFeature++;

                        if (currmatched)
                        {
                            goodMoras++;
                            break;
                        }
                    }
                }

                m_currScore = float(goodMoras) / float(m_currEntry.pitchcode.size());
                m_scoreText.setString(std::string("Accuracy: ") + std::to_string(int(m_currScore * 100)) + "%");
            }
        }

        auto calcPitchCircleRadius = [](float volume)
        {
            volume = std::max(volume, 0.1f);

            float volumeScale = volume;
            volumeScale /= 0.05f;
            if (volumeScale > 1.0f)
            {
                volumeScale = 1.0f;
            }

            return kPitchCircleRadius * volumeScale;
        };

        float volumeBeginTime = 0.0f;
        for (const auto& sample : m_currWordPitchHistory)
        {
            if (sample.volume > 0.01f)
            {
                volumeBeginTime = sample.time;
                break;
            }
        }

        auto calcPitchCircleX = [this, currWordTimeDuration, volumeBeginTime](float time)
        {
            float leftBound = kPitchCircleLeftBoundary + kInitialWordTime / kSecondsPerMora * kPitchCircleHorizontalSpacing;
            float x = kPitchCircleLeftBoundary + (time - volumeBeginTime) / (currWordTimeDuration - volumeBeginTime) * kPitchCircleHorizontalSpacing * float(m_currEntry.pitchcode.size());
            x = std::max(x, leftBound);
            return x;
        };

        float historyMinPitch = INFINITY;
        float historyMaxPitch = -INFINITY;
        if (m_currWordPitchHistory.empty())
        {
            historyMinPitch = kSpeechLowestFrequency;
            historyMaxPitch = kSpeechHighestFrequency;
        }
        else
        {
            for (const auto& sample : m_currWordPitchHistory)
            {
                historyMinPitch = std::min(historyMinPitch, sample.pitch);
                historyMaxPitch = std::max(historyMaxPitch, sample.pitch);
            }
        }

        auto calcPitchCircleY = [historyMinPitch, historyMaxPitch](float pitch)
        {
            float normalizedCurrPitch = (pitch - historyMinPitch) / (historyMaxPitch - historyMinPitch);
            if (std::isnan(normalizedCurrPitch))
                normalizedCurrPitch = 0.0f;
            if (normalizedCurrPitch < 0.0f)
                normalizedCurrPitch = 0.0f;
            if (normalizedCurrPitch > 1.0f)
                normalizedCurrPitch = 1.0f;
            return kPitchCircleTopBoundary - kPitchCircleHighToLowHeightDiff * normalizedCurrPitch;
        };

        m_pitchCircle.setOutlineThickness(kPitchCircleOutlineThickness);
        m_pitchCircle.setOutlineColor(sf::Color::Red);
        m_pitchCircle.setFillColor(sf::Color::Black);
        m_pitchCircle.setRadius(calcPitchCircleRadius(m_currVolume));
        m_pitchCircle.setOrigin(m_pitchCircle.getRadius(), m_pitchCircle.getRadius());
        m_pitchCircle.setPosition(calcPitchCircleX(m_currWordTime), calcPitchCircleY(m_currPitch));

        auto catmullRom = [](float p0, float p1, float p2, float p3, float t)
        {
            float t2 = t * t;
            float t3 = t2 * t;

            float x =
                0.5f * ((2.0f * p1) +
                (-p0 + p2) * t +
                (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);

            return x;
        };

        float lastSampleTime = m_currWordPitchHistory.empty() ? kInitialWordTime : m_currWordPitchHistory.back().time;
        const float kPitchPathVertsPerSecond = 100.0f;
        int numPitchPathVerts = (int)ceil((lastSampleTime - kInitialWordTime) * kPitchPathVertsPerSecond);
        std::vector<PitchHistorySample> pitchPathVerts(numPitchPathVerts);
        for (int i = 0; i < numPitchPathVerts; i++)
        {
            float t = kInitialWordTime + float(i) / kPitchPathVertsPerSecond;
            if (t <= kInitialWordTime) t = kInitialWordTime + 0.00001f;
            if (t >= lastSampleTime - kInitialWordTime) t = (lastSampleTime - kInitialWordTime) - 0.000001f;

            auto sampIts = std::equal_range(begin(m_currWordPitchHistory), end(m_currWordPitchHistory), PitchHistorySample{ t,0,0 },
                [](const PitchHistorySample& a, const PitchHistorySample& b) {
                return a.time < b.time;
            });
            assert(sampIts.second != end(m_currWordPitchHistory));

            sampIts.first = sampIts.first - 1;

            PitchHistorySample p0, p1, p2, p3;
            p1 = *sampIts.first;
            p2 = *sampIts.second;

            if (sampIts.first == begin(m_currWordPitchHistory))
            {
                p0.time = p1.time - (p2.time - p1.time);
                p0.pitch = p1.pitch - (p2.pitch - p1.pitch);
                p0.volume = p1.volume - (p2.volume - p1.volume);
            }
            else
            {
                p0 = *(sampIts.first - 1);
            }

            if (sampIts.second == end(m_currWordPitchHistory) - 1)
            {
                p3.time = p2.time + (p2.time - p1.time);
                p3.pitch = p2.pitch + (p2.pitch - p1.pitch);
                p3.volume = p2.volume + (p2.volume - p1.volume);
            }
            else
            {
                p3 = *(sampIts.second + 1);
            }

            PitchHistorySample interpolated;
            interpolated.time = t;

            float interpT;
            if (p2.time == p1.time)
            {
                interpT = 0.0f;
            }
            else
            {
                interpT = (t - p1.time) / (p2.time - p1.time);
            }

            interpolated.pitch = catmullRom(p0.pitch, p1.pitch, p2.pitch, p3.pitch, interpT);
            interpolated.volume = catmullRom(p0.volume, p1.volume, p2.volume, p3.volume, interpT);

            pitchPathVerts[i] = interpolated;
        }

        m_pitchCirclePathVerts.resize(numPitchPathVerts ? (numPitchPathVerts - 1) * 6 : 0);
        for (int i = 0; i < numPitchPathVerts - 1; i++)
        {
            float volume0 = pitchPathVerts[i].volume;
            float volume1 = pitchPathVerts[i + 1].volume;

            float time0 = pitchPathVerts[i].time;
            float time1 = pitchPathVerts[i + 1].time;

            float pitch0 = pitchPathVerts[i].pitch;
            float pitch1 = pitchPathVerts[i + 1].pitch;

            float radius0 = calcPitchCircleRadius(volume0);
            float left0 = calcPitchCircleX(time0) - radius0;
            float top0 = calcPitchCircleY(pitch0) - radius0;
            
            float radius1 = calcPitchCircleRadius(volume1);
            float left1 = calcPitchCircleX(time1) - radius1;
            float top1 = calcPitchCircleY(pitch1) - radius1;

            float x0 = left0 + radius0;
            float x1 = left1 + radius1;

            float y00 = top0 + radius0 * 2 * 0.8f;
            float y01 = top0 + radius0 * 2 * 0.2f;
            float y10 = top1 + radius1 * 2 * 0.8f;
            float y11 = top1 + radius1 * 2 * 0.2f;

            m_pitchCirclePathVerts[i * 6 + 0].position = sf::Vector2f(x0, y00);
            m_pitchCirclePathVerts[i * 6 + 1].position = sf::Vector2f(x0, y01);
            m_pitchCirclePathVerts[i * 6 + 2].position = sf::Vector2f(x1, y11);
            m_pitchCirclePathVerts[i * 6 + 3].position = sf::Vector2f(x0, y00);
            m_pitchCirclePathVerts[i * 6 + 4].position = sf::Vector2f(x1, y11);
            m_pitchCirclePathVerts[i * 6 + 5].position = sf::Vector2f(x1, y10);

            for (int j = 0; j < 6; j++)
            {
                m_pitchCirclePathVerts[i * 6 + j].color = sf::Color::Red;
            }
        }

        // Update the visualization of the pitch pattern
        {
            int numCircles = KatakanaOps::MoraCount(m_currEntry.katakana) + 1;
            
            bool odaka = false;
            if (!m_currEntry.pitchcode.empty() && (m_currEntry.pitchcode.back() & 0x8))
            {
                odaka = true;
            }

            m_guideCircles.resize(numCircles);
            for (int i = 0; i < numCircles; i++)
            {
                bool tailCircle = i == numCircles - 1;

                m_guideCircles[i].setOutlineThickness(kPitchCircleOutlineThickness);
                m_guideCircles[i].setOutlineColor(sf::Color::White);

                if (tailCircle)
                {
                    m_guideCircles[i].setFillColor(sf::Color::Black);
                }
                else
                {
                    m_guideCircles[i].setFillColor(sf::Color::White);
                }

                bool highAccent = false;
                if (tailCircle)
                {
                    if (odaka)
                    {
                        highAccent = false;
                    }
                    else if (m_currEntry.pitchcode.size() == 1)
                    {
                        highAccent = !(m_currEntry.pitchcode.back() & 0x1);
                    }
                    else
                    {
                        highAccent = m_currEntry.pitchcode.back() & 0x1;
                    }
                }
                else
                {
                    highAccent = m_currEntry.pitchcode[i] & 0x1;
                }

                m_guideCircles[i].setRadius(kPitchCircleRadius);
                m_guideCircles[i].setOrigin(kPitchCircleRadius, kPitchCircleRadius);
                m_guideCircles[i].setPosition(
                    kPitchCircleLeftBoundary + kPitchCircleHorizontalSpacing * i,
                    kPitchCircleTopBoundary - highAccent * kPitchCircleHighToLowHeightDiff);
            }

            m_guideCirclePathVerts.resize((m_guideCircles.size() - 1) * 6);
            for (int i = 0; i < numCircles - 1; i++)
            {
                float x0 = m_guideCircles[i].getGlobalBounds().left + m_guideCircles[i].getGlobalBounds().width / 2.0f;
                float x1 = m_guideCircles[i + 1].getGlobalBounds().left + m_guideCircles[i + 1].getGlobalBounds().width / 2.0f;

                float y00 = m_guideCircles[i].getGlobalBounds().top + m_guideCircles[i].getGlobalBounds().height * 0.6f;
                float y01 = m_guideCircles[i].getGlobalBounds().top + m_guideCircles[i].getGlobalBounds().height * 0.4f;
                float y10 = m_guideCircles[i + 1].getGlobalBounds().top + m_guideCircles[i + 1].getGlobalBounds().height * 0.6f;
                float y11 = m_guideCircles[i + 1].getGlobalBounds().top + m_guideCircles[i + 1].getGlobalBounds().height * 0.4f;

                m_guideCirclePathVerts[i * 6 + 0].position = sf::Vector2f(x0, y00);
                m_guideCirclePathVerts[i * 6 + 1].position = sf::Vector2f(x0, y01);
                m_guideCirclePathVerts[i * 6 + 2].position = sf::Vector2f(x1, y11);
                m_guideCirclePathVerts[i * 6 + 3].position = sf::Vector2f(x0, y00);
                m_guideCirclePathVerts[i * 6 + 4].position = sf::Vector2f(x1, y11);
                m_guideCirclePathVerts[i * 6 + 5].position = sf::Vector2f(x1, y10);

                for (int j = 0; j < 6; j++)
                {
                    m_guideCirclePathVerts[i * 6 + j].color = sf::Color::White;
                }
            }
        }

        // layout objects from bottom to top
        float currOffset = 0.0f;

        // Consider bottom padding
        currOffset += viewTopBottomPadding + (view.getSize().y - kPitchCircleHighToLowHeightDiff - kPitchCircleTopBoundary) - 130;

        // place kanji text
        float kanjiTextOffset = currOffset + m_kanjiText.getGlobalBounds().height + (m_kanjiText.getGlobalBounds().top - m_kanjiText.getPosition().y);
        currOffset = kanjiTextOffset;

        // Bit of padding between the texts
        currOffset += view.getSize().y * 0.02f;

        // place katakana text
        float worstKatakanaOffset = 0.0f;
        for (int i = 0; i < (int)m_katakanaTexts.size(); i++)
        {
            const sf::Text& text = m_katakanaTexts[i];
            float katakanaTextOffset = text.getGlobalBounds().height + (text.getGlobalBounds().top - text.getPosition().y);
            worstKatakanaOffset = std::max(worstKatakanaOffset, katakanaTextOffset);
        }
        for (int i = 0; i < (int)m_katakanaTexts.size(); i++)
        {
            sf::Text& text = m_katakanaTexts[i];
            text.setPosition(sf::Vector2f(
                kPitchCircleLeftBoundary + i * kPitchCircleHorizontalSpacing - text.getGlobalBounds().width / 2.0f, 
                view.getSize().y - currOffset - worstKatakanaOffset));
        }
        currOffset += worstKatakanaOffset;

        m_kanjiText.setPosition(sf::Vector2f(
            m_katakanaTexts[0].getGlobalBounds().left - (m_katakanaTexts[0].getGlobalBounds().left - m_katakanaTexts[0].getPosition().x),
            view.getSize().y - kanjiTextOffset));

        // Bit more paddings
        currOffset += view.getSize().y * 0.02f;

        // place romaji text
        float worstRomajiOffset = 0.0f;
        for (int i = 0; i < (int)m_romajiTexts.size(); i++)
        {
            const sf::Text& text = m_romajiTexts[i];
            float romajiTextOffset = text.getGlobalBounds().height + (text.getGlobalBounds().top - text.getPosition().y);
            worstRomajiOffset = std::max(romajiTextOffset, romajiTextOffset);
        }
        for (int i = 0; i < (int)m_romajiTexts.size(); i++)
        {
            sf::Text& text = m_romajiTexts[i];
            text.setPosition(sf::Vector2f(
                m_katakanaTexts[i].getGlobalBounds().left - (m_katakanaTexts[i].getGlobalBounds().left - m_katakanaTexts[i].getPosition().x),
                view.getSize().y - currOffset - worstRomajiOffset));
        }
        currOffset += worstRomajiOffset;
    }

    void draw(sf::RenderTarget& target, sf::RenderStates states) const override
    {
        if (!m_hasCurrEntry)
        {
            return;
        }

        target.draw(m_instructionsText);

        if (isnan(m_currScore))
        {
            target.draw(m_defaultScoreText);
        }
        else
        {
            target.draw(m_scoreText);
        }

        target.draw(m_kanjiText);
        
        for (const sf::Text& katakanaText : m_katakanaTexts)
        {
            target.draw(katakanaText);
        }
        for (const sf::Text& romajiText : m_romajiTexts)
        {
            target.draw(romajiText);
        }

        target.draw(m_guideCirclePathVerts);

        for (const auto& guideCircle : m_guideCircles)
        {
            target.draw(guideCircle);
        }

        target.draw(m_pitchCirclePathVerts);
        target.draw(m_pitchCircle);
    }

    void startStopGuideAnimation()
    {
        if (m_playGuideAnimation || (!m_playGuideAnimation && !m_currWordPitchHistory.empty()))
        {
            m_playGuideAnimation = false;
        }
        else
        {
            m_playGuideAnimation = true;
        }
        
        m_currWordTime = kInitialWordTime;
        m_currWordPitchHistory.clear();

        m_currScore = NAN;
    }
};

int main()
{
    // Create the main window
    sf::RenderWindow window(sf::VideoMode(1280, 720), "Pitch Accent Tutor");

    // Read the dictionary
    PitchAccentDictionary pitchDictionary;
    if (!pitchDictionary.loadFromFile("dictionary.csv"))
    {
        printf("Error: Failed to read dictionary.\n");
        return -1;
    }

    // Create a graphical text to display
    sf::Font font;
    if (!font.loadFromFile("meiryo.ttc"))
    {
        printf("Error: Failed to open font meiryo.ttc\n");
        return -1;
    }

    if (!VoiceRecorder::isAvailable())
    {
        printf("Error: Voice recording device not available.\n");
        return -1;
    }

    std::vector<std::string> availableDevices = VoiceRecorder::getAvailableDevices();
    printf("Voice recording device: %s\n", availableDevices.front().c_str());

    VoiceRecorder recorder;
    recorder.start();

    PitchGuideRenderer guideRenderer(font);
    
    int currDictionaryWord = 0;
    guideRenderer.setCurrentWord(pitchDictionary.at(currDictionaryWord));

    sf::Clock deltaClock;

    // Start the game loop
    while (window.isOpen())
    {
        sf::Time dt = deltaClock.restart();

        bool shouldUpdateCurrentWord = false;

        // Process events
        sf::Event event;
        while (window.pollEvent(event))
        {
            // Close window: exit
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::Left)
                {
                    currDictionaryWord = currDictionaryWord - 1;
                    if (currDictionaryWord == -1)
                    {
                        currDictionaryWord = (int)pitchDictionary.size() - 1;
                    }
                    shouldUpdateCurrentWord = true;
                }
                else if (event.key.code == sf::Keyboard::Right)
                {
                    currDictionaryWord = currDictionaryWord + 1;
                    currDictionaryWord = currDictionaryWord % (int)pitchDictionary.size();
                    shouldUpdateCurrentWord = true;
                }
                else if (event.key.code == sf::Keyboard::Space)
                {
                    guideRenderer.startStopGuideAnimation();
                }
            }
        }

        if (shouldUpdateCurrentWord)
        {
            guideRenderer.setCurrentWord(pitchDictionary.at(currDictionaryWord));
        }

        guideRenderer.setCurrentPitch(recorder.getCurrentPitch());
        guideRenderer.setCurrentVolume(recorder.getCurrentVolume());

        // Update the state of the program
        {
            guideRenderer.update(dt.asSeconds(), window);
        }

        // Render the current frame of animation
        {
            // Clear screen
            window.clear();

            // Draw the pitch accent guide
            window.draw(guideRenderer);

            // Update the window
            window.display();
        }
    }

    printf("Stopping recording...\n");

    recorder.stop();
}
