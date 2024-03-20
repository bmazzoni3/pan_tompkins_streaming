CC = gcc
CFLAGS = -O3 -g 

#Creare prerequisiti
pan_tompkins_streaming.x: pan_tompkins_streaming.o single_convolution.o
	$(CC) pan_tompkins_streaming.o single_convolution.o -o pan_tompkins_streaming.x  

#creare obiettivo: creazione file oggetto
pan_tompkins_streaming.o: pan_tompkins_streaming.c 
	$(CC) $(CFLAGS) -c pan_tompkins_streaming.c 

single_convolution.o: single_convolution.c 
	$(CC) $(CFLAGS) -c single_convolution.c 

#all è un altro obiettivo: creare file eseguibile
all: pan_tompkins_streaming.o single_convolution.o pan_tompkins_streaming.x 

clean: 
	rm -rf pan_tompkins_streaming.o single_convolution.o pan_tompkins_streaming.x || true

run:
	./pan_tompkins_streaming.x	


