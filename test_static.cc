#include <iostream>
#include <new>
#include "RankOptions.hpp"

char * nnew(const long row, const long measure_col_cnt, const long dim_col_cnt)
{
	long col = measure_col_cnt + dim_col_cnt;
	char *storage = new char[ row * col * sizeof(double) + row * sizeof(double) + 12]; 
	uint32_t* header = (uint32_t*) storage;
	header[0] = (uint32_t) row;
	header[1] = (uint32_t) measure_col_cnt;
	header[2] = (uint32_t) dim_col_cnt;
	header += 3;
	double** data = (double **)header;
	for (int i = 0; i < row; ++i )
	{
		data[i] = (double*)(storage + 3 * sizeof(uint32_t) + row * sizeof(double) + i * col * sizeof(double)); 
	};  
	return storage;
};



int main()
{
	const long row = 2000;
	const long mcol = 5;
	const long dcol = 5;
	char *data = nnew(row, mcol, dcol);
	uint32_t *header = (uint32_t*) data;
	std::cout << " row: " << header[0] << std::endl;
	std::cout << " mcol: " << header[1] << std::endl;
	std::cout << " dcol: " << header[2] << std::endl;
	header += 3;
	double ** real_data = (double**) header;
	srand(3);
	for (auto i = 0; i < row; ++i)
		for (auto j = 0; j < mcol + dcol; ++j)
			real_data[i][j] = rand() % 10000;

	std::vector<long> dim_cols = {5};
	std::vector<long> m_cols = {-1,-2};
	std::vector<unsigned long> r_cols = {5,1,2,3};
	std::vector<long> range_vecs;
	RankOptions rankOptions(dim_cols, m_cols, r_cols, 3, true, range_vecs, eRANK, data, 3, 1); 
	char * result = PSort::sort(rankOptions);
  
	header = (uint32_t*) result;
	std::cout << " result row: " << header[0] << std::endl;
	std::cout << " result mcol: " << header[1] << std::endl;
	std::cout << " result dcol: " << header[2] << std::endl;
	uint32_t result_row = header[0];
	uint32_t result_mcol = header[1];
	uint32_t result_dcol = header[2];
	header += 3;
	real_data = (double**) header;

	for (int i=0; i< result_row; ++i)
	{   
		for (int j=0; j< result_mcol + result_dcol; ++j)
		{   
			std::cout << real_data[i][j] << ",";
		}   
		std::cout << std::endl;
	};  
	delete [] result;

}
