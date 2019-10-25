CXX=g++
CXXFLAGS=-Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS+=-Wno-error=misleading-indentation
endif

ifndef HIJETMETHODS
$(error HIJETMETHODS is not set at all. Do 'source setJetMethodsEnv.sh'. Make terminating.)
endif

INCLUDE=-I$(HIJETMETHODS)
ROOT=`root-config --cflags --glibs`
FASTJET=`$(FASTJETPATH)/bin/fastjet-config --libs --cxxflags`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdf mkdirOutput bin/jetStatTesting.exe bin/plotJetStatTesting.exe bin/makeVacuumAreaHist.exe bin/evalPytHyd.exe bin/plotEvalPytHyd.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

bin/jetStatTesting.exe: src/jetStatTesting.C
	$(CXX) $(CXXFLAGS) src/jetStatTesting.C $(ROOT) $(FASTJET) $(INCLUDE) -o bin/jetStatTesting.exe

bin/plotJetStatTesting.exe: src/plotJetStatTesting.C
	$(CXX) $(CXXFLAGS) src/plotJetStatTesting.C $(ROOT) $(FASTJET) $(INCLUDE) -o bin/plotJetStatTesting.exe

bin/makeVacuumAreaHist.exe: src/makeVacuumAreaHist.C
	$(CXX) $(CXXFLAGS) src/makeVacuumAreaHist.C $(ROOT) $(FASTJET) $(INCLUDE) -o bin/makeVacuumAreaHist.exe

bin/evalPytHyd.exe: src/evalPytHyd.C
	$(CXX) $(CXXFLAGS) src/evalPytHyd.C $(ROOT) $(FASTJET) $(INCLUDE) -o bin/evalPytHyd.exe

bin/plotEvalPytHyd.exe: src/plotEvalPytHyd.C
	$(CXX) $(CXXFLAGS) src/plotEvalPytHyd.C $(ROOT) $(FASTJET) $(INCLUDE) -o bin/plotEvalPytHyd.exe

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
	rm -f src/*~
	rm -f src/#*#
	rm -f src/*.so
	rm -f src/*.pcm
	rm -f src/*.d
