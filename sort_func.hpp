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
enum SORT_OPERATOR 
{TB_COUNT=0,TB_SUM,TB_PERCENT,RANK,DENSERANK,NTILE};

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

class Util
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
		static float ** mark_count_rank(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result = nullptr)
		{
			int idx = index; double current_value=0; int current_type = 1;
			while (idx <= size)
			{
				double tmp = func_(data[idx], formula_args_);
				if (tmp != current_value)
				{
					current_value= tmp;
					current_type = idx - index + 1;
				}
				data[idx][result_idx_] = current_type;
				idx++;
			};
			return data;
		};
		static float **  mark_denserank(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result = nullptr)
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

			return data;
		};
		static float **  mark_ntile(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result = nullptr)
		{
			//top_bottom_parms: [int1, int2]
			int part = (size - index + 1) / top_bottom_param[0];
			int remains = (size - index + 1) % top_bottom_param[0];
			int idx = index; double current_value=0; int current_type = 0;
			while (idx <= size)
			{
				current_value ++;
				if (current_type == top_bottom_param[1])
				{
					data[idx][result_idx_] = current_type;
					//memcpy(result[0], data[idx], col_size_ * sizeof(float));
					//safer when memory overlaps
					memmove(result[0], data[idx], col_size_ * sizeof(float));
					result++;
				};
				if (current_value == part + (remains > 0))
				{
					//data[idx][result_idx_] = current_type;
					current_value = 0;
					current_type++;
					if ( current_type > top_bottom_param[1] ) 
					{
						return result;
					};
					remains--;
				};
				idx++;
			};
			return result;
		};
		// so far, just one sum value considered
		static float **  mark_sum(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result = nullptr)
		{
			//top_bottom_param: [sum1, sum2] -- array of sums
			int idx = index; double running_total=0; int current_type = 0;
			while (idx <= size)
			{
				running_total += func_(data[idx], formula_args_);
				data[idx][result_idx_] = current_type;
				//memcpy(result[0], data[idx], col_size_ * sizeof(float));
				memmove(result[0], data[idx], col_size_ * sizeof(float));
				result++;				
				if ( running_total > top_bottom_param[0] )
				{
					return result;
				};
				idx++;
			};
			return result;
		};
		static float **  mark_percent(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result = nullptr)
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
				return data;
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
			return data;
		};
		typedef float ** (*OP_FUNC)(float **data, int size, int index, double sum, const int result_idx_, const int *formula_args_, float(*func_)(const float*data, const int*args), const int* top_bottom_param, const int top_bottom_psize_, int col_size_, float ** result);

		static void top_bottom_sum(float** data, const int data_size_, std::vector<double>& sum, std::vector<int>& index,const int *formula_args_, float (*func_)(const float* data, const int* args), const int dim_sort_size_, const int* sort_, const SORT_OPERATOR tb_type )
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

		static float ** new_new(const int row, const int col)
		{
			float ** data = new float* [row];
			for (int i = 0; i < row; ++i)
				data[i] = new float[col];
			return data;
		};
		static int new_delete(float** data, const int row, const int col, float** end)
		{
			int final_data_size = end - data;
			for (int i = final_data_size; i < row; ++i)
				delete [] data[i];
			return final_data_size;
		};

	public:
		static constexpr OP_FUNC operator_mapping[6]  = {mark_count_rank, mark_sum, mark_percent, mark_count_rank, mark_denserank, mark_ntile};
		// parallel sorting function, std::sort under certain threashold
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
			static int sort(float** data, const int data_size_, const int size_, const int measure_idx_, const int dim_sort_size_,
					const int *sort_, const int *formula_args_, float (*func_)(const float* data, const int* args),
					const int result_idx_, const bool asc, const int* top_bottom_param, const int top_bottom_psize_, S sort_type, const SORT_OPERATOR tb_type, int max_thread=0)
			{
				assert(data != NULL);
				assert(data_size_ != 0);
				assert(size_ != 0);
				assert(measure_idx_ != 0);
				assert(func_ != NULL);
				assert(result_idx_ >0);
				// calling sort
				Functor cmp(size_, measure_idx_, dim_sort_size_, sort_, formula_args_, func_, asc);
				sort_type(data, data+data_size_, cmp, max_thread);
				// summing for each partition
				std::vector<double> sum; std::vector<int> index;
				top_bottom_sum(data, data_size_, sum, index, formula_args_, func_, dim_sort_size_, sort_, tb_type);

				// filling results
				int idx = 0;
				float ** inter_data = data;
				for(int i = 0; i < index.size(); ++i)
				{
					inter_data = Util::operator_mapping[tb_type](data, index[i], idx, sum[i], result_idx_, formula_args_, func_, top_bottom_param, top_bottom_psize_, size_, inter_data); 
					idx = index[i] + 1;
				}

				if ( tb_type == NTILE || tb_type == TB_SUM )
				{
					return new_delete(data, data_size_, size_, inter_data);
				}else
				{
					return data_size_;
				};
			};
};


/*********end of function/class **********/
constexpr Util::OP_FUNC Util::operator_mapping[6];

