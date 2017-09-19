all: sub_main_finder tree_filesplitter_plussubs find_orphans tree_filesplitter_orphans sub_main_finder_MassCut_ParticleNum_orphans
orphans: tree_filesplitter_orphans sub_main_finder_MassCut_ParticleNum_orphans find_orphans
main: sub_main_finder tree_filesplitter_plussubs
fosters: tree_filesplitter_fosters

sub_main_finder:
	gcc -o sub_main_finder sub_main_finder.c -Wall -O3

sub_main_finder_MassCut_ParticleNum_orphans:
	gcc -o sub_main_finder_MassCut_ParticleNum_orphans sub_main_finder_MassCut_ParticleNum_orphans.c -Wall -O3

tree_filesplitter_plussubs:
	gcc -o tree_filesplitter_plussubs tree_filesplitter_plussubs.c -Wall -O3

tree_filesplitter_orphans: 
	gcc -o tree_filesplitter_orphans tree_filesplitter_orphans.c -Wall -O3

tree_filesplitter_fosters:
	gcc -o tree_filesplitter_fosters tree_filesplitter_fosters.c -Wall -O3

find_orphans:
	gcc -o find_orphans find_orphans.c -Wall -O3

.PHONY: clean

clean :
	rm -f sub_main_finder
	rm -f tree_filesplitter_plussubs
	rm -f sub_main_finder_MassCut_ParticleNum_orphans
	rm -f tree_filesplitter_orphans
	rm -f find_orphans
	rm tree_filesplitter_fosters
