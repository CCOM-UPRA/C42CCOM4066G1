# Included Libraries

cxxopts.hpp is a library utilized to handle the input of  parameters, and their storage, to the program.


This library was selected from the following based on the given criteria. Libraries were found by ChatGPT, and further researched for the final selection.
| **Library** | **Pros** | **Cons** |
|-------------|----------|---------|
| [**cxxopts**](https://github.com/jarro2783/cxxopts) | - Very clean and simple<br>- Good defaults<br>- Excellent documentation and examples<br>- Recently maintained and up to date<br>- Grouped help messages and custom help | - Limited support for deeply nested subcommands<br>- Slight learning curve for advanced features |
| [**argparse (C++)**](https://github.com/p-ranav/argparse) | - Header-only, easy integration<br>- Familiar API inspired by Python’s argparse<br>- Supports sub-commands and default values<br>- Actively maintained with frequent updates | - Requires C++17 or newer<br>- Smaller community and fewer resources<br>- No built-in config or env parsing<br>- Slightly increases compile times for large projects |
| [**CLI11**](https://github.com/CLIUtils/CLI11) | - Very rich feature set<br>- Header-only, C++11 compatible<br>- Actively developed and well-documented<br>- Widely adopted in the C++ community | - Larger header, slightly longer compile times<br>- Overkill for simple programs |
| [**TCLAP**](https://github.com/idokka/tclap) | - Header-only and lightweight<br>- Mature and stable<br>- Works with older C++ standards<br>- Detailed manual and examples available | - Largely unmaintained<br>- No subcommand support<br>- Non-standard parsing behaviors<br>- Less convenient value retrieval |

Github Repository (Needed File): [cxxopts.hpp](https://github.com/jarro2783/cxxopts/blob/master/include/cxxopts.hpp)

# Additional Information

###  Removed Stackable Short Flag
The stackable short flag was removed because in c++ (using our current library) it caused compatibility issues with the SAM short flag, given they were both -s.

### Variable Name Changes
Not implemented for now, some variable names currently have typos (such as 'treshold'), and they dont seem to follow a standard...

### Running command (for test cases)
```
g++ -std=c++17  "Aidan G. Varela Soto/precision_recall.cpp" -o precision_recall 
./precision_recall -g genome.fa -f probes.fasta -s sam_dir -c meta_dir -l 1000
```

