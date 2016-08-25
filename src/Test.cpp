#include "cute.h"
#include "ide_listener.h"
#include "xml_listener.h"
#include "cute_runner.h"
#include "kdtree.h"
#include <set>
#include <array>
#include <cstdio>
#include "kdtree.cpp"

using namespace dhlib;
using namespace std;

template class KDTree<3,int>;
template class KDTree<4,string>;

void testEmptyTree(){
	KDTree<3,int> tree; // a 3d tree holding ints
	ASSERTM("tree is empty so begin should equal end", tree.begin() == tree.end());
}

void testRootLeaf(){
	array<string, 8> strings = {"one","two","three","four","five","six","seven","eight"};
	KDTree<4,string> tree;
	tree(1,2,3,4) = strings[0];
	tree(5,6,7,8) = strings[1];
	tree(0,0,0,0) = strings[2];
	tree(1,-1,2,3) = strings[3];
	tree(6,-3,3,6) = strings[4];
	tree(11,2,3,88,0) = strings[5];
	tree(9,4,2,6) = strings[6];
	tree(1,1,1,1) = strings[7];
	set<string> elems (tree.values.begin(), tree.values.end());
	int count = 0;
	for(auto str = strings.begin(); str != strings.end(); ++str){
		char message[256];
		sprintf(message, "tree contains string %s", str->c_str());
		ASSERTM(message, elems.find(*str) != elems.end());
		count++;
	}
	ASSERTM("tree has 8 elements", count == 8);
}

bool runAllTests(int argc, char const *argv[]) {
	cute::suite s { };
	s.push_back(CUTE(testEmptyTree));
	s.push_back(CUTE(testRootLeaf));
	cute::xml_file_opener xmlfile(argc, argv);
	cute::xml_listener<cute::ide_listener<>> lis(xmlfile.out);
	auto runner { cute::makeRunner(lis, argc, argv) };
	bool success = runner(s, "AllTests");
	return success;
}

int main(int argc, char const *argv[]) {
    return runAllTests(argc, argv) ? EXIT_SUCCESS : EXIT_FAILURE;
}
