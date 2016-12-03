# Utilisation :
#
# - Compiler le code:
#
#     $ make
#
# - Executer le code:
#
#     $ make run
#
# - Supprimer les sorties (binaires, données)
#
#     $ make clean
#

.PHONY: all run clean

all: main

main: main.cpp
	g++ -o $@ $<

temperature.gnu: main
	./$<

run: temperature.gnu

clean:
	rm -f main temperature.gnu

