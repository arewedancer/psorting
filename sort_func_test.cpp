#include "sort_func.hpp"
#include "gtest/gtest.h"
#include <algorithm>
#include <iostream>
using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestCase;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;


#define data_size 11000000
#define col_size 5

class TestFixture : public ::testing::Test
{
	protected:
	float ** data;
	int row_size; 
	TestFixture()
	{
		row_size = data_size;
	};
	void SetUp()
	{
		srand(3);
		data = new float*[data_size];
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

		std::cout << "first 10 lines of test data:" << std::endl;
		std::for_each(data, data+10, [](float* e)
				{
				std::for_each(e, e+col_size-1, [](float d){
						std::cout << std::fixed << std::setprecision(3) << d << ",";});
			  std::cout << std::endl;	
				});

	};
	void TearDown()
	{
		std::cout << "after sorting, first 10 lines of test data:" << std::endl;
		std::for_each(data, data+10, [](float* e)
				{
				std::for_each(e, e+col_size, [](float d){
						std::cout << std::fixed << std::setprecision(3) << d << ",";});
			  std::cout << std::endl;	
				});
		for(int i=0; i< row_size; ++i)
			delete [] data[i];
		delete [] data;
	};

	~TestFixture()
	{
		/*for(int i=0; i< data_size; ++i)
			delete [] data[i];
		delete [] data;*/

	};

};


TEST_F(TestFixture, TBB_MERGE_Rank)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {30, 40, 30};
	int return_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_tbb<float*, Functor>, RANK);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};
TEST_F(TestFixture, Cilk_MERGE_Rank)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {30, 40, 30};
	int return_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_cilkplus<float*, Functor>, RANK);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};
TEST_F(TestFixture, SAMPLE_Rank)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {30, 40, 30};
	int return_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_sample_sort<float*, Functor>, RANK);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};

TEST_F(TestFixture, TBB_MERGE_TB_PERCENT)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {50,50};
	int return_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_tbb<float*, Functor>, TB_PERCENT);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};

TEST_F(TestFixture, TBB_MERGE_NTILE)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {3,1};
	row_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_tbb<float*, Functor>, NTILE);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};

TEST_F(TestFixture, TBB_MERGE_SUM)
{
	tbb::tick_count t0 = tbb::tick_count::now();
	int args[]={3};
	int sort_size = 1;
	int sort[] = {0};
	int top_bottom_args[] = {30000};
	row_size = Util::sort(data, data_size, 5, 3, sort_size, sort, args, measure_func, 4, false, top_bottom_args, 3, parallel_merge_sort_tbb<float*, Functor>, TB_SUM);
	tbb::tick_count t1 = tbb::tick_count::now();
  std::cout << "time taken: " << (t1-t0).seconds() << std::endl;
};

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


