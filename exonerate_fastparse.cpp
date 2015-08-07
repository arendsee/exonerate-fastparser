/* exonerate_fastparse.cpp
 * 
 * Reads raw exonerate output, currently only parsing the vulgar line
 *
 */

#define BUFFER_SIZE 4096

#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
using namespace std;

void parse_vulgar(string line, bool stop);

int main()
{
    char buffer[BUFFER_SIZE];
    string line;
    int roll = 0;
    bool between = true;
    bool stop = false;
    while(fgets(buffer, BUFFER_SIZE, stdin)){
        line = buffer;
        if(between){
            if(line.length() > 3 && line.substr(2,6) == "Target"){
                roll = 1; 
                between = false;
            }
            continue;
        }

        switch(roll % 5) {
            case 2: if(line[0] == 'v'){
                        parse_vulgar(line, stop);
                        between = true;
                        stop = false;
                    }
            case 4: if(line.find('*') != size_t(-1))
                        stop = true;
        }
        roll++;
    }
}

void parse_vulgar(string line, bool stop){
    stringstream s(line);
    string word;
    for(int i = 0; i < 10 && s >> word; i++){
        if(i > 0){
            cout << word << '\t';
        }
    }
    cout << stop << endl;
    for(int i = 0; s >> word; i++){
        if(i % 3 == 0)
            cout << "> ";
        cout << word;
        if(i % 3 != 2){
            cout << '\t';
        } else {
            cout << endl;
        }
    }
}
