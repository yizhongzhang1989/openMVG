//
// Created by v-xinli1 on 10/28/2020.
//

#include <fstream>
#include <vector>

int main(int argc, char **argv)
{
    std::string input_file(argv[1]);
    std::string output_file(argv[2]);


    std::ifstream fin(input_file, std::ifstream::in);
    if( !fin.is_open() ) return -1;
    std::ofstream fout( output_file, std::ofstream::out );
    fout.precision(20);
    int i=0;
    while( !fin.eof() )
    {
        double time;
        fin >> time;
//        if( i++%2 == 0 )
        {
            double new_time;
            new_time = 1403715273262142976 + time*1e9;
            fout << new_time << std::endl;
        }
    }
    fin.close();
    fout.close();


    return 0;
}

//1403715273262142976
//1403715273
//262142976