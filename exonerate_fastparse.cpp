/* exonerate_fastparse.cpp
 * 
 * Reads raw exonerate output, currently only parsing the vulgar line
 *
 */

#define BUFFER_SIZE 4096

#include <cstdlib>
#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
using namespace std;

void parse_vulgar_simple( string line, int stop, char delim);
void parse_vulgar_fulltab(string line, int stop, char delim);
void parse_vulgar_verbose(string line, int stop, char delim);
int get_first_integer(string line);
int get_stop_position(string line, int offset);
void write_fulltab_header();

int main()
{
    write_fulltab_header();

    // buffer to hold data input
    char buffer[BUFFER_SIZE];

    // string to hold line input
    string line;

    // index for ordering within a record
    int position = 0;

    // flag for being within the alignment region of a record
    bool between = true;

    // query start position for current alignment segment
    int running_position = 0;

    // position of first stop codon, 0 indicates none found
    int stop = 0;

    // temporary storage for output of get_stop_position
    int tmp_stop;

    while(fgets(buffer, BUFFER_SIZE, stdin)){

        // cast the cstring as a string object
        line = buffer;

        if(between){
            if(line.length() > 3 && line.substr(2,6) == "Target"){
                position = 1; 
                between = false;
            }
            continue;
        }

        switch(position % 5) {
            case 2: if(line[0] == 'v'){
                        parse_vulgar_fulltab(line, stop, '\t');
                        between = true;
                        stop = 0;
                    } else {
                        running_position = get_first_integer(line);
                    }
            case 4: if(stop == 0) {
                        tmp_stop = get_stop_position(line, running_position);
                        if(tmp_stop != 0 && stop == 0){
                            stop = tmp_stop; 
                        }
                    }
        }

        position++;
    }
}

int get_stop_position(string line, int pos){
    for(unsigned int i = 0; i < line.length(); i++){
        // All capital letters signal a new amino acid
        if(line[i] >= 'A' && line[i] <= 'Z'){
            pos++;            
        }
        // I am interested only in the position of the first stop,
        // so return when first '*' is encountered
        else if(line[i] == '*'){
            return(pos);
        }
    }
    // Return 0 if no stop is found
    return 0;
}

int get_first_integer(string line){
    int x;
    std::stringstream ss;
    ss << line;
    ss >> x;
    return x;
}

void write_fulltab_header(){
    cout << "query_seqid\t"
         << "query_start\t"
         << "query_stop\t"
         << "query_strand\t"
         << "target_contig\t"
         << "target_start\t"
         << "target_stop\t"
         << "target_strand\t"
         << "score\t"
         << "first_stop\t"
         << "has_frameshift\t"
         << "num_split_codons\t"
         << "num_intron\t"
         << "max_intron"
         << endl;
}

void parse_vulgar_fulltab(string line, int stop, char delim){
    stringstream s(line);
    string word;
    for(int i = 0; i < 10 && s >> word; i++){
        if(i > 0){
            cout << word << delim;
        }
    }
    cout << stop << delim;
    bool shifted = false;
    int nintrons = 0;
    int splits = 0;
    int longest_intron = 0;
    while(s >> word){
        switch(word[0]){
            case 'F': shifted = true;
            case 'I': { nintrons++;
                        s >> word;
                        s >> word;
                        int len = atoi(word.c_str());
                        if(len > longest_intron)
                            longest_intron = len;
                      }
            case 'S': splits++;
        }
    }
    cout << shifted        << delim
         << splits         << delim
         << nintrons       << delim
         << longest_intron << endl;
}

void parse_vulgar_simple(string line, bool stop, char delim){
    stringstream s(line);
    string word;
    for(int i = 0; i < 9 && s >> word; i++){
        if(i > 0){
            cout << word << delim;
        }
    }
    s >> word;
    cout << word << endl;
}

void parse_vulgar_verbose(string line, bool stop, char delim){
    stringstream s(line);
    string word;
    for(int i = 0; i < 10 && s >> word; i++){
        if(i > 0){
            cout << word << delim;
        }
    }
    cout << stop << endl;
    for(int i = 0; s >> word; i++){
        if(i % 3 == 0)
            cout << "> ";
        cout << word;
        if(i % 3 != 2){
            cout << delim;
        } else {
            cout << endl;
        }
    }
}
