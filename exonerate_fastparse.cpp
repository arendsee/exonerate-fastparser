/* exonerate_fastparse.cpp
 * 
 * Reads raw exonerate output (from STDIN)
 *
 */

#define BUFFER_SIZE 4096

#include <cstdlib>
#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <getopt.h>
using namespace std;

string get_sans_introns_header();
string get_with_introns_header();
string get_simple_header();

string parse_simple(string line, int stop, char delim);
string parse_sans_introns(string line, int stop, char delim);
string parse_with_introns(string line, int stop, char delim);
string parse_verbose(string line, int stop, char delim);
int get_first_integer(string line);
int get_stop_position(string line, int offset);

int main(int argc, char **argv)
{

    int format = 0;
    int c;

    while (1) {
        static struct option long_options[] =
        {
            {"simple",       no_argument, &format, 0}, 
            {"sans-introns", no_argument, &format, 1}, 
            {"with-introns", no_argument, &format, 2}, 
            {"verbose",      no_argument, &format, 3},
            {"help",         no_argument, 0, 'h'},
            {0,0,0,0}
        };
        int option_index = 0;
        c = getopt_long (argc, argv, "h", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch(c) {
            case 0:
                break;
            case 'h':
                cout << "Usage: exonerate-fastparse [options] filename" << endl
                     << "Format options:" << endl
                     << "   --simple - writes tabular output with the following columns:" << endl
                     << "      1. query_seqid"   << endl
                     << "      2. query_start"   << endl
                     << "      3. query_stop"    << endl
                     << "      4. query_strand"  << endl
                     << "      5. target_seqid"  << endl
                     << "      6. target_start"  << endl
                     << "      7. target_stop"   << endl
                     << "      8. target_strand" << endl
                     << "      9. score"         << endl
                     << "   --sans-introns - adds 5 columns to the simple output" << endl
                     << "      10. first_stop"     << endl
                     << "      11. has_frameshift" << endl
                     << "      12. split_codons"   << endl
                     << "      13. introns"        << endl
                     << "      14. max_introns"    << endl
                     << "   --with-introns - adds 1 column to the sans-introns output" << endl
                     << "      15. intron_lengths - a comma-delimited list of intron lengths" << endl
                     << "   --verbose\tprints simple tabular output and then one line for each gff entry" << endl;
                return 0;
        }
    }

    string (*parse)(string, int, char);

    switch(format) {
        case 0:
            cout << get_simple_header() << endl;
            parse = parse_simple;
            break;
        case 1:
            cout << get_sans_introns_header() << endl;
            parse = parse_sans_introns;
            break;
        case 2:
            cout << get_with_introns_header() << endl;
            parse = parse_with_introns;
            break;
        case 3:
            parse = parse_verbose;
            break;
    }

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
                        cout << parse(line, stop, '\t') << endl;
                        between = true;
                        stop = 0;
                    } else {
                        running_position = get_first_integer(line);
                    }
                    break;
            case 4: if(stop == 0) {
                        tmp_stop = get_stop_position(line, running_position);
                        if(tmp_stop != 0 && stop == 0){
                            stop = tmp_stop; 
                        }
                    }
                    break;
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

string get_simple_header(){
    string header = "query_seqid\t"      
                    "query_start\t"      
                    "query_stop\t"       
                    "query_strand\t"     
                    "target_contig\t"    
                    "target_start\t"     
                    "target_stop\t"      
                    "target_strand\t"    
                    "score\t";
    return header;
}

string get_sans_introns_header(){
    string header = get_simple_header();
    header += "first_stop\t"       
              "has_frameshift\t"   
              "split_codons\t" 
              "introns\t"       
              "max_intron";
    return header;
}

string get_with_introns_header(){
    string header = get_sans_introns_header();
    header += "\tintron_lengths";
    return header;
}

string parse_sans_introns(string line, int stop, char delim){
    stringstream s(line);
    stringstream out;
    string word;
    for(int i = 0; i < 10 && s >> word; i++){
        if(i > 0){
            out << word << delim;
        }
    }
    out << stop << delim;
    bool shifted = false;
    int nintrons = 0;
    int splits = 0;
    int longest_intron = 0;
    while(s >> word){
        switch(word[0]){
            case 'F': shifted = true;
                      break;
            case 'I': { nintrons++;
                        s >> word;
                        s >> word;
                        int len = atoi(word.c_str());
                        if(len > longest_intron)
                            longest_intron = len;
                      }
                      break;
            case 'S': splits++;
                      break;
        }
    }
    out << shifted  << delim
        << splits   << delim
        << nintrons << delim
        << longest_intron;
    return out.str();
}

string parse_with_introns(string line, int stop, char delim){
    stringstream out;
    stringstream s(line);
    out << parse_sans_introns(line, stop, delim) << delim;
    string word;
    bool has_intron = false;
    while(s >> word){
        if(word[0] == 'I'){
           s >> word;
           s >> word;
           if(has_intron)
               out << ',';
           out << word;
           has_intron = true; 
        }
    }
    if(! has_intron)
        out << '-';
    return out.str();
}

string parse_simple(string line, int stop, char delim){
    stringstream s(line);
    stringstream out;
    string word;
    for(int i = 0; i < 9 && s >> word; i++){
        if(i > 0){
            out << word << delim;
        }
    }
    s >> word;
    out << word;
    return out.str();
}

string parse_verbose(string line, int stop, char delim){
    stringstream s(line);
    stringstream out;
    string word;
    for(int i = 0; i < 10 && s >> word; i++){
        if(i > 0){
            out << word << delim;
        }
    }
    out << stop << endl;
    for(int i = 0; s >> word; i++){
        if(i % 3 == 0)
            out << "> ";
        out << word;
        out << ((i % 3 != 2) ? delim  : '\n');
    }
    return out.str();
}
