

To run:

1) Type "make" to create all necessary executables. GCC is the default compiler used, but the Makefile can easily be edited to use a different compiler. Python codes were written in 2.7, but should run fine up to version 3.6 (at least), although it is reccommend to run the test file first to confirm outputs are correct.

2) Run the code using: ./get_halo_props.sh path_to_tree_files tree_#_#_# tree_#_#_#

	- path_to_tree_files is the path (relative or absolute) to where the tree_#_#_#.dat files can be accessed. The files in this directory should be ORDERED, by being run through Manodeep's Sinha's publicly available code sort_forests (https://github.com/manodeep/sort_forests) first. Otherwise, some of the data points may be missing.

	- tree_#_#_# are the range of files to be searched. (Ex tree_0_0_0 tree_0_0_9 as inputs would search and find subhalos in the first 10 files, tree_0_0_0, tree_0_0_1, tree_0_0_2,.... up to and including tree_0_0_9)

	- output will be one file, called subhalo_properties.txt that will be created in the directory in which the code was run


TO RUN THE TEST FILE (reccommended if using a different version of python than 2.7, if the tree files may have been modified at any point, or just if you want to quickly confirm output before running code for 2 days) :


Some notes:

	- Processing for the trees will take some time: it is highly reccommended to open multiple screen sessions on a computer and run the data is smaller chunks. Expect the processing of ~20 trees to take between 5-10 hours.

	- Tree files are drastically different sizes, some will take up to ~5x as long to process as others.
