CXX=g++
CXXFLAGS=-Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS+=-Wno-error=misleading-indentation
endif

ifndef HIJETMETHODS
$(error HIJETMETHODS is not set at all. Do 'source setJetMethodsEnv.sh'. Make terminating.)
endif

INCLUDE=-I$(HIJETMETHODS)
LIB=-L$(HIJETMETHODS)/MakePytHyd/lib
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_LIB=mkdir -p $(PWD)/lib
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirLib mkdirPdf mkdirOutput lib/libconfigParser.so lib/libgeneralTreeHandler.so bin/combinePythiaWithBkgd.exe bin/testConfigParser.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

lib/libconfigParser.so: src/configParser.C
	$(CXX) $(CXXFLAGS) -shared -fPIC src/configParser.C $(ROOT) $(INCLUDE) -o lib/libconfigParser.so

lib/libgeneralTreeHandler.so: src/generalTreeHandler.C
	$(CXX) $(CXXFLAGS) -shared -fPIC src/generalTreeHandler.C $(ROOT) -lEG $(INCLUDE) $(LIB) -lconfigParser -o lib/libgeneralTreeHandler.so

bin/combinePythiaWithBkgd.exe: src/combinePythiaWithBkgd.C
	$(CXX) $(CXXFLAGS) src/combinePythiaWithBkgd.C $(ROOT) -lEG $(INCLUDE) $(LIB) -lconfigParser -lgeneralTreeHandler -o bin/combinePythiaWithBkgd.exe

bin/testConfigParser.exe: src/testConfigParser.C
	$(CXX) $(CXXFLAGS) src/testConfigParser.C $(ROOT) -lEG $(INCLUDE) $(LIB) -lconfigParser -lgeneralTreeHandler -o bin/testConfigParser.exe

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f config/*~
	rm -f config/#*#
	rm -f include/*~
	rm -f include/#*#
	rm -f input/*~
	rm -f input/#*#
	rm -f lib/*.so
	rm -rf lib
	rm -f src/*~
	rm -f src/#*#
	rm -f src/*.so
	rm -f src/*.pcm
	rm -f src/*.d
