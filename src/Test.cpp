#include "cute.h"
#include "ide_listener.h"
#include "xml_listener.h"
#include "cute_runner.h"
#include "kdtree.h"
#include "kdtree.cpp"

using namespace dhlib;

template class KDTree<3,int>;

void testEmptyTree(){
	KDTree<3,int> tree; // a 3d tree holding ints
	ASSERTM("tree is empty so begin should equal end", tree.begin() == tree.end());
}

void testRootLeaf(){
	KDTree<4,string> tree;
}

bool runAllTests(int argc, char const *argv[]) {
	cute::suite s { };
	s.push_back(CUTE(testEmptyTree));
	cute::xml_file_opener xmlfile(argc, argv);
	cute::xml_listener<cute::ide_listener<>> lis(xmlfile.out);
	auto runner { cute::makeRunner(lis, argc, argv) };
	bool success = runner(s, "AllTests");
	return success;
}

int main(int argc, char const *argv[]) {
    return runAllTests(argc, argv) ? EXIT_SUCCESS : EXIT_FAILURE;
}
