#include <iostream>
#include <algorithm>
#include <iomanip>
#include <unistd.h>
#include <cassert>
#include <vector>
#include <memory>
#include <cmath>
#include "tbb/tbb.h"
#include "merge_sort_tbb.h"
#include "merge_tbb.h"
//#include "compat/thread"
#include <thread>
#include "sample_sort_util.h"
#include "quicksort_util.h"
#include "quicksort_tbb_task.h"
#include "map_keys_to_bins_tbb.h"
#include "bin_tbb.h"
#include "repack_and_subsort_tbb.h"
#include "sample_sort.h"
#include <cilk/cilk.h>
#include "serial_merge.h"
#include "merge_cilk.h"
#include "merge_sort_cilk.h"
#define data_size 11000000 
#define col_size 5
//#define col_size 1
#define max_thread 100
enum TOP_BOTTOM_TYPE
{COUNT=0,SUM,PERCENT,RANK,DENSERANK,NTILE};

// measure function: return measure value
static float measure_func(const float* data, const int* args)
{
	//std::cout << "func " << args[0] << std::endl;
	return data[args[0]];
};

// compare function object for sorting
struct Functor {
	int size;
	int dim_sort_size; 
	int measure_sort_size; 
	int measure_idx; 
	const int* sort;
	const int* formula_args;
	float (*func)(const float*, const int* args);
	bool asc;
	explicit Functor(const int size_, const int measure_idx_, const int dim_sort_size_, 
			const int *sort_, const int *formula_args_, float (*func_)(const float*, const int* args), const bool asc_):
		sort(sort_), measure_idx(measure_idx_), formula_args(formula_args_), size(size_),
		dim_sort_size(dim_sort_size_), func(func_),asc(asc_)
	{};  
	bool operator() (const float* a, const float* b) const { 
		for (int i=0; i< dim_sort_size; ++i )
		{
			if ( a[sort[i]] == b[sort[i]]) continue;
			return (a[sort[i]] > b[sort[i]]) ^ asc;
		}
		return (func(a, formula_args) > func(b, formula_args)) ^ asc;

	};
};

class util
{
	private:
		static bool is_same_partition(const float *a, const float*b, const int sort_size, const int* sort)
		{
			for(int i=0; i<sort_size; ++i)
			{
				if (a[sort[i]] != b[sort[i]])  return false;
			}
			return true;
		};

//refactoring: replace switch-case with array of function pointers
		static void mark_count_rank(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_)
		{
						int idx = index; double current_value=0; int current_type = 1;
						while (idx <= size)
						{
							double tmp = func_(data[idx], formula_args_);
							if (tmp != current_value)
							{
								current_value= tmp;
								current_type = idx + 1;
							}
							data[idx][result_idx_] = current_type;
							idx++;
						};
		};
		static void mark_denserank(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_)
		{
						int idx = index; double current_value=0; int current_type = 0;
						while (idx <= size)
						{
							double tmp = func_(data[idx], formula_args_);
							if (tmp != current_value)
							{
								current_value= tmp;
								current_type++;
							}
							data[idx][result_idx_] = current_type;
							idx++;
						};
		};
		static void mark_ntile(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_)
		{
						//top_bottom_parms: [int1, int2]
						int part = (size - index + 1) / top_bottom_param[0];
						int remains = (size - index + 1) % top_bottom_param[0];
						int idx = index; double current_value=0; int current_type = 0;
						while (idx <= size)
						{
							current_value ++;
							if (current_value == part + (remains > 0))
							{
								data[idx][result_idx_] = current_type;
								current_value = 0;
								current_type++;
								remains--;
							}else {
								data[idx][result_idx_] = current_type;
							};
							idx++;
						};
		};
		static void mark_sum(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_)
		{
						int idx = index; double running_total=0; int current_type = 0;
						while (idx <= size)
						{
							if (current_type == -1)
							{
								data[idx][result_idx_] = current_type;
								idx++;
								continue;
							};
							running_total += func_(data[idx], formula_args_);
							if ( running_total > top_bottom_param[current_type] )
							{
								if (++current_type <= top_bottom_psize_ -1 )
								{
									running_total = func_(data[idx], formula_args_);
								}else{
									current_type = -1;
								};
							};

							data[idx][result_idx_] = current_type;
							idx++;
						};
		};
		static void mark_percent(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_)
		{
						//corner case size - index + 1 < top_bottom_psize_
						if ( size - index + 1 <= top_bottom_psize_)
						{
							int idx = index;
							while (idx <= size)
							{
								data[idx][result_idx_] = 0;
								idx++;

							};
							return;
						}

						int idx = index; double running_total=0; int current_type = 0;
						while (idx <= size)
						{
							running_total += func_(data[idx], formula_args_);
							if ( running_total > top_bottom_param[current_type] * sum/100 && current_type != top_bottom_psize_-1 )
							{
								running_total = func_(data[idx], formula_args_);
								current_type++;
							};

							// mark
							data[idx][result_idx_] = current_type;
							idx++;
						};
		};
		typedef void (*OP_FUNC)(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_);
		//static constexpr OP_FUNC operator_mapping[6]  = {mark_count_rank, mark_sum, mark_percent, mark_count_rank, mark_denserank, mark_ntile};

		static void mark_partition(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, const TOP_BOTTOM_TYPE tb_type)
		{
			//TODO: validate parameters
			switch(tb_type)
			{
				//for tie-breaker, rank is the same, but count++
				//so A,B,C, A, B both 1, C 3
				case COUNT:
				case RANK:
					{
						int idx = index; double current_value=0; int current_type = 1;
						while (idx <= size)
						{
							double tmp = func_(data[idx], formula_args_);
							if (tmp != current_value)
							{
								current_value= tmp;
								current_type = idx + 1;
							}
							data[idx][result_idx_] = current_type;
							idx++;
						};
						break;
					}
				case DENSERANK:
					{
						int idx = index; double current_value=0; int current_type = 0;
						while (idx <= size)
						{
							double tmp = func_(data[idx], formula_args_);
							if (tmp != current_value)
							{
								current_value= tmp;
								current_type++;
							}
							data[idx][result_idx_] = current_type;
							idx++;
						};
						break;
					}

				case NTILE:
					{
						//top_bottom_parms: [int1, int2]
						int part = (size - index + 1) / top_bottom_param[0];
						int remains = (size - index + 1) % top_bottom_param[0];
						int idx = index; double current_value=0; int current_type = 0;
						while (idx <= size)
						{
							current_value ++;
							if (current_value == part + (remains > 0))
							{
								data[idx][result_idx_] = current_type;
								current_value = 0;
								current_type++;
								remains--;
							}else {
								data[idx][result_idx_] = current_type;
							};
							idx++;
						};
						break;
					}
				case SUM:
					{
						int idx = index; double running_total=0; int current_type = 0;
						while (idx <= size)
						{
							if (current_type == -1)
							{
								data[idx][result_idx_] = current_type;
								idx++;
								continue;
							};
							running_total += func_(data[idx], formula_args_);
							if ( running_total > top_bottom_param[current_type] )
							{
								if (++current_type <= top_bottom_psize_ -1 )
								{
									running_total = func_(data[idx], formula_args_);
								}else{
									current_type = -1;
								};
							};

							data[idx][result_idx_] = current_type;
							idx++;
						};
						break;
					}
				case PERCENT:
					{
						//corner case size - index + 1 < top_bottom_psize_
						if ( size - index + 1 <= top_bottom_psize_)
						{
							int idx = index;
							while (idx <= size)
							{
								data[idx][result_idx_] = 0;
								idx++;

							};
							return;
						}

						int idx = index; double running_total=0; int current_type = 0;
						while (idx <= size)
						{
							running_total += func_(data[idx], formula_args_);
							if ( running_total > top_bottom_param[current_type] * sum/100 && current_type != top_bottom_psize_-1 )
							{
								running_total = func_(data[idx], formula_args_);
								current_type++;
							};

							// mark
							data[idx][result_idx_] = current_type;
							idx++;
						};
					}
			}

		};

		static void top_bottom_sum(float** data, const int data_size_, std::vector<double>& sum, std::vector<int>& index,const int *formula_args_, float (*func_)(const float* data, const int* args), const int dim_sort_size_, const int* sort_, const TOP_BOTTOM_TYPE tb_type )
		{
			int idx = 1; double tmp_sum = func_(data[0], formula_args_);
			while (idx < data_size_)
			{
				if (!is_same_partition(data[idx], data[idx-1], dim_sort_size_, sort_))
				{
					index.emplace_back(idx-1);
					sum.emplace_back(tmp_sum);
					tmp_sum = 0;
				};
				tmp_sum += func_(data[idx], formula_args_);
				idx++;
			}
			index.emplace_back(idx-1);sum.emplace_back(tmp_sum);
		};

	public:
		//typedef void (*OP_FUNC)(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_);
		static constexpr OP_FUNC operator_mapping[6]  = {mark_count_rank, mark_sum, mark_percent, mark_count_rank, mark_denserank, mark_ntile};
		// TOP/BOTTOM function
		// function parameters:
		// S sort_type [parallel_merge_sort_tbb, parallel_merge_sort_cilkplus, parallel_sample_sort]
		// tb_type: operators of [COUNT, SUM, PERCENT]
		// float ** data: data to be operated on
		// data_size_: size of rows
		// size_: size of columns
		// measure_idx_: index of measure in a row
		// dim_sort_size_: partition column size
		// sort_: pointer to partition column index
		// func_: measure operations
		// formula_args_: measure operation argument pointer
		// result_idx_: result column index
		// asc: true- BOTTOM, false- TOP
		// top_bottom_param: top_bottom_parameter pointer, eg. {20,30,50}
		// top_bottom_psize_: top_botton_parameter size
		template<typename S>
			static void top_bottom(float** data, const int data_size_, const int size_, const int measure_idx_, const int dim_sort_size_,
					const int *sort_, const int *formula_args_, float (*func_)(const float* data, const int* args),
					const int result_idx_, const bool asc, const int* top_bottom_param, const int top_bottom_psize_, S sort_type, const TOP_BOTTOM_TYPE tb_type)
			{
				assert(data != NULL);
				assert(data_size_ != 0);
				assert(size_ != 0);
				assert(measure_idx_ != 0);
				assert(func_ != NULL);
				assert(result_idx_ >0);
				// calling sort
				Functor cmp(size_, measure_idx_, dim_sort_size_, sort_, formula_args_, func_, asc);
				sort_type(data, data+data_size_, cmp);
				// summing for each partition
				std::vector<double> sum; std::vector<int> index;
				top_bottom_sum(data, data_size_, sum, index, formula_args_, func_, dim_sort_size_, sort_, tb_type);

				// filling results
				int idx = 0;
				for(int i = 0; i < index.size(); ++i)
				{
					//mark_partition(data, index[i], idx, sum[i], result_idx_, formula_args_, func_, top_bottom_param, top_bottom_psize_, tb_type); 
					util::operator_mapping[tb_type](data, index[i], idx, sum[i], result_idx_, formula_args_, func_, top_bottom_param, top_bottom_psize_); 
					idx = index[i] + 1;
				}
			};


};
constexpr util::OP_FUNC util::operator_mapping[6];

void thread_func(float ** data)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	//int top_bottom_args[] = {30, 40, 30};
	//int top_bottom_args[] = {300, 400, 300};
	int top_bottom_args[] = {4,1};
	util::top_bottom(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_tbb<float*, Functor>, NTILE); 
	tbb::tick_count t1 = tbb::tick_count::now();
	std::cout << (t1-t0).seconds() << std::endl;
};

//#define data_size 11000000
//#define data_size 1100
//#define col_size 5
//#define max_thread 100
int main()
{
	srand(3);
	float ** data = new float*[data_size];
	//std::unique_ptr<float*> data(new float*[data_size]);
	for (int i=0; i< data_size; ++i)
	{
		float * row = new float[col_size]; 
		//std::unique_ptr<float> row(new float[col_size]); 
		//for (int j=0; j<col_size-1; ++j)
		for (int j=0; j<col_size-1; ++j)
		{
			//row[j] = rand() % 50000;
			row[j] = rand() % 100000;
			//std::cout << "data " << row[j] << std::endl;
		}
		data[i] = row;
	};
	/*std::thread threads[max_thread];
		/int args[]={3};
		int sort_size = 1;
		int sort[] = {0};
		Functor cmp(5, 3, sort_size, sort, args, measure_func, true);
	//std::cout << parallel_sample_sort_sampling(data, data+data_size, cmp, true) << std::endl;
	//parallel_sample_sort_sampling(data, data+data_size, cmp, true);
	//std::sort(data, data+data_size, cmp);
	parallel_sample_sort<float*, Functor>(data, data+data_size, cmp, true);
	exit(1);
	for ( int i=0; i< max_thread; ++i)
	{
	//T* test = new T[size];
	//    int tsize = random ? rand()%size : size;
	float ** new_data = new float*[data_size];
	std::copy(data, data+data_size, new_data); 
	threads[i] = std::thread(thread_func, new_data);
	};  
	for ( int i=0; i< max_thread; ++i)
	{
	threads[i].join();
	};*/ 


	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	//int top_bottom_args[] = {30, 40, 30};
	//int top_bottom_args[] = {300, 400, 300};
	int top_bottom_args[] = {4,1};
	//Less less(5, 3, 0, NULL, args, measure_func);

	for (int i=0; i< 10; ++i)
	{
		for (int j=0; j<col_size; ++j)
			std::cout << std::fixed << std::setprecision(5) << data[i][j] << ",";
		std::cout << std::endl;
	};
	std::cout << "sorting" << std::endl;
	thread_func(data);

	for (int i=0; i< 2000; ++i)
	{
		for (int j=0; j<col_size; ++j)
			std::cout << std::fixed << std::setprecision(5) << data[i][j] << ",";
		std::cout << std::endl;
	};

};
