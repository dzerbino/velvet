README.TXT

VELVET SOURCE 
March 28 2008
Daniel Zerbino

NOTE: The PDF manual in this directory contains all the information contained
in this text file, plus much more!

> SUMMARY
	* A/ REQUIREMENTS
	* B/ COMPILING INSTRUCTIONS

----------------------------------------------------------------------------------
A/ REQUIREMENTS

	Velvet should function on any standard 64bit Linx environment with
gcc. A good amount of physical memory (12GB to start with, more is no luxury)
is recommended. 

----------------------------------------------------------------------------------
B/ COMPILING INSTRUCTIONS

Normally, with a GNU environment, just type:

> make

For colorspace Velvet replace that command with 

> make color

Otherwise compile each *.c file separately, then execute the default
instructions at the top of Makefile. 
