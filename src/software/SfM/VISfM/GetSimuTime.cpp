//
// Created by v-xinli1 on 10/28/2020.
//

#include <fstream>
#include <vector>

int main(int argc, char **argv)
{
    /*std::string input_file(argv[1]);
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
        if( i++%2 == 0 )
        {
            fout << time << std::endl;
        }
    }
    fin.close();
    fout.close();*/

    std::string output_file(argv[1]);
    std::ofstream fout( output_file, std::ofstream::out );
    fout.precision(3);
    double step = 1. / 20.;
    for( int i=0;i<200;++i )
    {
        fout << static_cast<double>( i ) * step << std::endl;
    }
    fout.close();

    return 0;
}