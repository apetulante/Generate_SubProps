all: sub_main_finder tree_filesplitter_plussubs full_split_halos

sub_main_finder:
	gcc -o sub_main_finder sub_main_finder.c -O3

tree_filesplitter_plussubs:
	gcc -o tree_filesplitter_plussubs tree_filesplitter_plussubs.c -Wall -O3

full_split_halos: 
	gcc -o full_split_halos full_split_halos.c -Wall -O3

.PHONY: clean

clean :
	rm -f sub_main_finder
	rm -f tree_filesplitter_plussubs
	rm -f full_split_halos
