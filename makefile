all: kmergo

KMERGO_BIN_DIR = bin
KMC_API_DIR = kmc_api
KMERGO_DIR = kmergo
STATS_ML_DIR = stats_ml

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11
CLINK	= -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11

STATS_ML_OBJS = \
$(STATS_ML_DIR)/stats.o \
$(STATS_ML_DIR)/logistic_classifier.o

KMERGO_OBJS = \
$(KMERGO_DIR)/kmer_counting.o \
$(KMERGO_DIR)/kmer_union.o \
$(KMERGO_DIR)/kmer_filtering.o \
$(KMERGO_DIR)/kmer_assembly.o \
$(KMERGO_DIR)/csv_lib.o \
$(KMERGO_DIR)/kmer_losertree.o \
$(KMERGO_DIR)/KmerGO.o

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

$(KMERGO_OBJS) $(KMC_API_OBJS) $(STATS_ML_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

kmergo: $(KMERGO_OBJS) $(KMC_API_OBJS) $(STATS_ML_OBJS)
	-mkdir -p $(KMERGO_BIN_DIR)
	$(CC) $(CLINK) -o $(KMERGO_BIN_DIR)/$@ $^

clean:
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(KMERGO_DIR)/*.o
	-rm -f $(STATS_ML_OBJS)/*.o
	-rm -rf bin
