CXX = mpicxx
CXXFLAGS = -isystem /opt/impi-5.1.3.223/intel64/include
CXXFLAGS += -O0 -mfpmath=sse -fstack-protector-all -g
CXXFLAGS += -W -Wall -Wextra -Wunused -Wcast-align -Werror
CXXFLAGS += -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith
CXXFLAGS += -Wformat-security -Wmissing-format-attribute -Wformat=1
CXXFLAGS += -Wwrite-strings -Wcast-align -Wno-long-long
CXXFLAGS += -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual
CXXFLAGS += -Wno-suggest-attribute=format -Wno-error=cast-function-type

SOURCES = Amethods.cpp Discrepancy.cpp Multiplication.cpp Rotation.cpp \
          main.cpp my_main.cpp number_translation.cpp Time.cpp
HEADERS = header.h

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = a.out

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $^ -o $@

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
