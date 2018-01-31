# JP Pitch Accent Tutor

## How to build and run

1. Download and install [Visual Studio](https://www.visualstudio.com/downloads/) (The free community edition is OK)
2. In the Visual Studio installer, check "Desktop development with C++" to get the C++ compiler.
3. Open `jp-pitch-accent-tutor.sln`
4. **Important**: Switch to [Release build](https://i.imgur.com/4HkIZNM.png)
5. Right click on the `jp-pitch-accent-tutor` project in the Solution Explorer, and choose ["Set as StartUp Project"](https://i.imgur.com/8mv8vOD.png).
6. Press "Local Windows Debugger" (the green play arrow) or press F5.

## How to add words to the dictionary

1. Go to `jp-pitch-accent-tutor/dictionary-builder`
2. Open `words_to_export.csv` in a text editor.
3. Add more lines to the file.

Notes:

* Each line must follow the format: `<kanji><TAB><katakana>`
	* OR the following format: `<kanji><TAB><katakana><TAB><disambiguation_index>`
* There must only be one TAB. Not spaces, not more than one tab.
* *Beware of text editors that automatically insert spaces when you press tab!!*
* The kanji and katakana must match *exactly* the corresponding entry in the nhk dictionary.
* The disambiguation index is for when the dictionary has more than one entry for the same word.
	* The disambiguation index says which one it should use.
	* Note the disambiguation index is 0-based, so the first entry in the dictionary is 0, the second is 1, etc.

The dictionary used by the pitch accent tutor will be automatically updated when you build it **in Release mode**. That's why the build steps say that it's so important to build in Release mode (at least once.)

