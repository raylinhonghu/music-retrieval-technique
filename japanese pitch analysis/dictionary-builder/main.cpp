#include <fcntl.h>
#include <io.h>

#include <cstdio>
#include <fstream>
#include <codecvt>

#include <unordered_set>
#include <unordered_map>
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

// See above.
std::wofstream OpenOutputUTF8FileAsWide(const char* filename)
{
    std::wofstream wof(filename);
    wof.imbue(std::locale(std::locale::empty(), new std::codecvt_utf8<wchar_t>));
    return std::move(wof);
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        printf("Usage: dictionary-builder <nhk_pronunciation.csv> <words_to_export.csv> <output.csv>\n");
        return -1;
    }

    const char* nhk_pronunciation_filename = argv[1];
    const char* words_to_export_filename = argv[2];
    const char* output_filename = argv[3];

    // Allow printing UTF-16 characters to the console.
    _setmode(_fileno(stdout), _O_U16TEXT);
    _setmode(_fileno(stderr), _O_U16TEXT);

    // The whole NHK pronunciation dictionary.
    std::wifstream nhk_pronunciation = OpenInputUTF8FileAsWide(nhk_pronunciation_filename);
    if (!nhk_pronunciation)
    {
        fprintf(stderr, "Error: Failed to open %s\n", nhk_pronunciation_filename);
        return -1;
    }

    // The words we care about to export.
    std::wifstream words_to_export = OpenInputUTF8FileAsWide(words_to_export_filename);
    if (!words_to_export)
    {
        fprintf(stderr, "Error: Failed to open %s\n", words_to_export_filename);
        return -1;
    }

    // The output file.
    std::wofstream output = OpenOutputUTF8FileAsWide(output_filename);
    if (!output)
    {
        fprintf(stderr, "Error: Failed to open %s\n", output_filename);
        return -1;
    }

    // Read all the words we care to export into a hash table.
    std::unordered_set<std::wstring> export_word_set;
    std::unordered_map<std::wstring, int> export_index;
    std::unordered_map<std::wstring, std::pair<int, int>> export_disambiguation;
    {
        std::wstring curr_word;
        int curr_word_index = 0;

        std::wstring line;
        while (std::getline(words_to_export, line))
        {
            if (line.find(' ') != std::wstring::npos)
            {
                fwprintf(stderr, L"Error: \"%s\" has space in it (should use tabs)\n", line.c_str());
                return -1;
            }

            size_t curr_token_start = 0, curr_token_end = 0;

            size_t kanji_katakana_start = curr_token_end;

            // Find the token for the kanji of the word
            curr_token_end = line.find('\t', curr_token_start);
            if (curr_token_end == std::wstring::npos)
            {
                continue;
            }
            curr_token_start = curr_token_end + 1;

            // Find the token for the katakana of the word
            curr_token_end = line.find('\t', curr_token_start);
            if (curr_token_end == std::wstring::npos)
            {
                curr_token_end = line.size();
            }
            // Read out the word's kanji and katakana
            curr_word.assign(line.c_str() + kanji_katakana_start, curr_token_end - kanji_katakana_start);
            curr_token_start = curr_token_end + 1;

            bool ambiguous = !export_word_set.insert(curr_word).second;
            if (ambiguous)
            {
                fwprintf(stderr, L"Error: Duplicated word %s\n", curr_word.c_str());
                return -1;
            }

            export_index.emplace(curr_word, curr_word_index);
            curr_word_index += 1;

            // Get the optional disambiguation number
            if (curr_token_end < line.size())
            {
                curr_token_end = line.find('\t', curr_token_start);
                if (curr_token_end == std::wstring::npos)
                {
                    curr_token_end = line.size();
                }
                int disambiguation = std::stoi(line.substr(curr_token_start, curr_token_end - curr_token_start));
                export_disambiguation.emplace(curr_word, std::make_pair(0, disambiguation)).second;
                curr_token_start = curr_token_end + 1;
            }
        }
    }

    std::unordered_set<std::wstring> found_word_set;

    std::vector<std::wstring> lines_to_output(export_index.size());

    // Find the pronunciation for all the dictionary words we care about
    {
        std::wstring curr_word;
        std::wstring curr_pitch_code;
        
        int curr_line_number = 0;

        std::wstring line;
        while (std::getline(nhk_pronunciation, line))
        {
            curr_line_number += 1;

            size_t curr_token_start = 0, curr_token_end = 0;

            // Extract the word at the start of this dictionary entry.
            {
                size_t kanji_katakana_start = curr_token_end;

                // Find the token for the kanji of the word
                curr_token_end = line.find('\t', curr_token_start);
                if (curr_token_end == std::wstring::npos)
                {
                    continue;
                }
                curr_token_start = curr_token_end + 1;

                // Find the token for the katakana of the word
                curr_token_end = line.find('\t', curr_token_start);
                if (curr_token_end == std::wstring::npos)
                {
                    continue;
                }
                // Read out the word's kanji and katakana
                curr_word.assign(line.c_str() + kanji_katakana_start, curr_token_end - kanji_katakana_start);
                curr_token_start = curr_token_end + 1;
            }

            // If this isn't a word we care about, skip it.
            if (export_word_set.find(curr_word) == end(export_word_set))
            {
                continue;
            }

            auto found_disambiguation = export_disambiguation.find(curr_word);
           
            if (found_disambiguation != end(export_disambiguation))
            {
                found_disambiguation->second.first += 1;

                if (found_disambiguation->second.first - 1 != found_disambiguation->second.second)
                {
                    // this wasn't the disambiguation we were looking for
                    continue;
                }
            }

            bool duplicate_word = !found_word_set.insert(curr_word).second;
            if (duplicate_word)
            {
                fwprintf(stderr, L"Error: Duplicate word \"%s\".\n", curr_word.c_str());
                return -1;
            }

            // Read the pitch accent info and encode the pitch code
            {
                curr_pitch_code.clear();

                // encoder state machine state
                bool currently_overline = false;
                bool currently_nasal = false;
                bool currently_nopron = false;

                bool* span_stack[3] = { 0,0,0 };
                int span_stack_sz = 0;

                // state IDs
                enum
                {
                    state_reading_chars,
                    state_starting_span,
                    state_reading_class_name,
                    state_finished_reading_class_name,
                    state_ending_span,
                    state_reading_unicode,
                };

                int starting_span_progress = 0;
                int reading_class_name_progress = 0;
                int reading_class_name_bits = 0x7;
                int ending_span_progress = 0;
                int reading_unicode_progress = 0;
                int reading_unicode_bits = 0x3;

                bool is_odaka = false;

                int curr_state = state_reading_chars;

                const char* parser_error = 0;

                for (size_t i = curr_token_start, sz = line.size(); i < sz; i++)
                {
                    wchar_t ch = line[i];

                    if (ch == '<')
                    {
                        if (curr_state != state_reading_chars)
                        {
                            parser_error = "invalid state transition due to <";
                            goto end_parse;
                        }

                        curr_state = state_starting_span;
                        starting_span_progress = 0;
                    }
                    else if (ch == '>')
                    {
                        if (!(curr_state == state_ending_span && ending_span_progress == 4) &&
                            !(curr_state == state_finished_reading_class_name))
                        {
                            parser_error = "unexpected >";
                            goto end_parse;
                        }

                        if (curr_state == state_ending_span)
                        {
                            span_stack_sz -= 1;
                            *span_stack[span_stack_sz] = false;
                        }

                        curr_state = state_reading_chars;
                    }
                    else if (ch == '/')
                    {
                        if (curr_state != state_starting_span)
                        {
                            parser_error = "invalid state transition due to /";
                            goto end_parse;
                        }

                        if (span_stack_sz == 0)
                        {
                            parser_error = "missing opening span";
                            goto end_parse;
                        }

                        curr_state = state_ending_span;
                        ending_span_progress = 0;
                    }
                    else if (ch == '&')
                    {
                        if (curr_state != state_reading_chars)
                        {
                            parser_error = "invalid state transition due to &";
                            goto end_parse;
                        }

                        curr_state = state_reading_unicode;
                        reading_unicode_progress = 0;
                        reading_unicode_bits = 0x3;
                    }
                    else if (curr_state == state_starting_span)
                    {
                        static const char starting_span_class[] = "span class=\"";
                        if (starting_span_progress < sizeof(starting_span_class) - 1)
                        {
                            if (ch == starting_span_class[starting_span_progress])
                            {
                                starting_span_progress += 1;
                                if (starting_span_progress == sizeof(starting_span_class) - 1)
                                {
                                    curr_state = state_reading_class_name;
                                    reading_class_name_progress = 0;
                                    reading_class_name_bits = 0x7;
                                }
                            }
                            else
                            {
                                parser_error = "expected span class";
                                goto end_parse;
                            }
                        }
                    }
                    else if (curr_state == state_reading_class_name)
                    {
                        static const char overline_class[] = "overline\"";
                        static const char nasal_class[] = "nasal\"";
                        static const char nopron_class[] = "nopron\"";

                        if (reading_class_name_bits & 0x1)
                        {
                            if (reading_class_name_progress >= sizeof(overline_class) - 1 ||
                                ch != overline_class[reading_class_name_progress])
                            {
                                reading_class_name_bits &= ~0x1;
                            }
                        }
                        if (reading_class_name_bits & 0x2)
                        {
                            if (reading_class_name_progress >= sizeof(nasal_class) - 1 ||
                                ch != nasal_class[reading_class_name_progress])
                            {
                                reading_class_name_bits &= ~0x2;
                            }
                        }
                        if (reading_class_name_bits & 0x4)
                        {
                            if (reading_class_name_progress >= sizeof(nopron_class) - 1 ||
                                ch != nopron_class[reading_class_name_progress])
                            {
                                reading_class_name_bits &= ~0x4;
                            }
                        }

                        reading_class_name_progress += 1;

                        if (reading_class_name_bits == 0)
                        {
                            parser_error = "unknown span class name";
                            goto end_parse;
                        }

                        bool* p_activated_flag = 0;
                        if (reading_class_name_bits == 0x1 && reading_class_name_progress == sizeof(overline_class) - 1)
                        {
                            p_activated_flag = &currently_overline;
                        }
                        else if (reading_class_name_bits == 0x2 && reading_class_name_progress == sizeof(nasal_class) - 1)
                        {
                            p_activated_flag = &currently_nasal;
                        }
                        else if (reading_class_name_bits == 0x4 && reading_class_name_progress == sizeof(nopron_class) - 1)
                        {
                            p_activated_flag = &currently_nopron;
                        }

                        if (p_activated_flag)
                        {
                            if (*p_activated_flag)
                            {
                                parser_error = "identical nested span class";
                                goto end_parse;
                            }

                            assert(span_stack_sz <= sizeof(span_stack) / sizeof(*span_stack));
                            span_stack[span_stack_sz] = p_activated_flag;
                            span_stack_sz += 1;

                            *p_activated_flag = true;

                            curr_state = state_finished_reading_class_name;
                        }
                    }
                    else if (curr_state == state_ending_span)
                    {
                        static const char span[] = "span";

                        if (ending_span_progress < sizeof(span) - 1 && ch == span[ending_span_progress])
                        {
                            ending_span_progress += 1;
                        }
                        else
                        {
                            parser_error = "expected </span>";
                            goto end_parse;
                        }
                    }
                    else if (curr_state == state_reading_unicode)
                    {
                        static const char downshift_accent[] = "#42780;";
                        static const char nasal_degree[] = "#176;";

                        if (reading_unicode_bits & 0x1)
                        {
                            if (reading_unicode_progress >= sizeof(downshift_accent) - 1 ||
                                ch != downshift_accent[reading_unicode_progress])
                            {
                                reading_unicode_bits &= ~0x1;
                            }
                        }
                        if (reading_unicode_bits & 0x2)
                        {
                            if (reading_unicode_progress >= sizeof(nasal_degree) - 1 ||
                                ch != nasal_degree[reading_unicode_progress])
                            {
                                reading_unicode_bits &= ~0x2;
                            }
                        }

                        reading_unicode_progress += 1;

                        if (reading_unicode_bits == 0)
                        {
                            parser_error = "unknown unicode bits";
                            goto end_parse;
                        }

                        bool is_end = true;

                        if (reading_unicode_bits == 0x1 && reading_unicode_progress == sizeof(downshift_accent) - 1)
                        {
                            if (i + 1 == sz)
                            {
                                is_odaka = true;
                            }
                        }
                        else if (reading_unicode_bits == 0x2 && reading_unicode_progress == sizeof(nasal_degree) - 1)
                        {
                            // do nothing
                        }
                        else
                        {
                            is_end = false;
                        }

                        if (is_end)
                        {
                            curr_state = state_reading_chars;
                        }
                    }
                    else if (curr_state == state_reading_chars)
                    {
                        int pitch_code_bits = 0;
                        if (currently_overline) pitch_code_bits |= 0x1;
                        if (currently_nasal) pitch_code_bits |= 0x2;
                        if (currently_nopron) pitch_code_bits |= 0x4;

                        curr_pitch_code.push_back(pitch_code_bits[L"0123456789ABCDEF"]);
                    }
                }

                if (curr_state != state_reading_chars)
                {
                    parser_error = "unexpected end of line";
                    goto end_parse;
                }

                if (is_odaka)
                {
                    curr_pitch_code.back() |= 0x8;
                }

                end_parse:
                
                if (parser_error)
                {
                    fprintf(stderr, "Error: Parser at line %d: %s\n", curr_line_number, parser_error);
                    return -1;
                }
            }

            // Output the kanji/katakana/pitchcode of the word
            lines_to_output[export_index.at(curr_word)] = curr_word + L"\t" + curr_pitch_code + L"\n";
        }
    }

    if (export_word_set.size() != found_word_set.size())
    {
        for (const auto& export_word : export_word_set)
        {
            if (found_word_set.find(export_word) == end(found_word_set))
            {
                fwprintf(stderr, L"Warning: Did not find word \"%s\" in dictionary.\n", export_word.c_str());
                return -1;
            }
        }
    }

    for (const std::wstring& line : lines_to_output)
    {
        output << line;
    }
}
