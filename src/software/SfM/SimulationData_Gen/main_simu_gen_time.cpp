//
// Created by v-xinli1 on 11/1/2020.
//

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <io.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

class CsvReader_Simu
{
public:
    CsvReader_Simu() = delete;
    CsvReader_Simu( const std::string& filename, const char split = ' ' )
    {
        std::cout << "open file" << filename << std::endl;
        split_= split;
        csvInput_.open( filename );
        if( !csvInput_ ) throw std::runtime_error( "Error csv file dict: " + filename );

        std::string header;
        getline( csvInput_, header );

        data_ = std::vector<double>(8, 0);
        last_data_ = std::vector<double>(8, 0);
    }
    ~CsvReader_Simu()
    {
        if(csvInput_) csvInput_.close();
    }

    bool readline()
    {
        if( !csvInput_ ) throw std::runtime_error(" Not Set File to Read ");
        if( csvInput_.eof() ) return false;

        std::string line;
        getline(csvInput_, line);

        std::istringstream readStr( line );
        std::string part;

        for( int i= 0;i<8;++i  )
        {
            getline( readStr, part, split_ );
            data_[i] = std::strtod( part.c_str(), nullptr );
        }

        data_[0] *= 1e3;
        data_[0] = static_cast<long long int>(data_[0]);

        if( last_data_[0] != 0 )
        {
            if( (data_[0] - last_data_[0]) != 5 )
            {
                std::cout << line << std::endl;
                std::cout << data_[0] << " - " <<last_data_[0] << " != 5" << std::endl;
                return false;
            }
        }
        last_data_ = data_;
        return true;
    }

    std::vector<double> data_;
private:
    std::vector<double> last_data_;
    std::ifstream csvInput_;
    char split_;
};


int main(int argc, char **argv)
{
    std::string input_file(argv[1]);
    std::string output_file(argv[2]);


    std::ofstream fout( output_file, std::ofstream::out );
    CsvReader_Simu reader( input_file );
    std::vector<double> times;
    while(reader.readline())
    {
        times.push_back(reader.data_[0]);
    }


    fout.precision(30);
    for( int i=0;i<times.size();++i )
    {
        if(i%20 == 0)
        {
            fout << times[i]*1e6 << std::endl;
        }
    }


    return 0;
}