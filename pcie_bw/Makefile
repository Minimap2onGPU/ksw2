SOURCES = test.cpp 
OBJECTS = $(SOURCES:.cpp=.o) 
EXECUTABLE = test

CXX=hipcc
LDFLAGS=
CXXFLAGS= -std=c++17 -Wno-unused-result -fopenmp -O2 #-O0 -g #--save-temps

all: ${EXECUTABLE} 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXECUTABLE): $(OBJECTS) 	
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDFLAGS) 

clean: 	
	rm -f $(EXECUTABLE) 	
	rm -f $(OBJECTS)
	rm -f *amdgcn* *x86_64* results*

