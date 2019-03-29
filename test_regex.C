#include <iostream>
#include <string>
#include <regex>

using namespace std;

int main() {
    string input = "position Si 33.1621 46.0566 -199.094";
    string output = std::regex_replace(
        input,
        regex("[^0-9]*([0-9]+).*"),
        string("\\1")
        );
    
    cout << input << endl;
    cout << output << endl;

}
