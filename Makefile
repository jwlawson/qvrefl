NAME = qvrefl
MAJOR = 0
MINOR = 1
VERSION = $(MAJOR).$(MINOR)

ifeq ($(CXX),icpc)
CXXFLAGS += -Wall -xHOST
OPT += -O3 -ipo
B_OPT += $(OPT)
AR = xiar
else
CXXFLAGS += -Wall -Wextra -march=native
CXXFLAGS += -fno-signed-zeros -fno-math-errno -fno-rounding-math
CXXFLAGS += -fno-signaling-nans -fno-trapping-math
CXXFLAGS += -ffinite-math-only -Wno-misleading-indentation
OPT += -Ofast
endif
CXXFLAGS += -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -DNDEBUG

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')

ifeq ($(uname_O),Cygwin)
	CXXFLAGS += -std=gnu++14 -DCYGWIN_STOI
endif
ifeq ($(uname_S),Linux)
	CXXFLAGS += -std=c++14
endif

TEST = test$(NAME)

BASE_DIR = .
SRC_DIR = $(BASE_DIR)/src
TEST_DIR = $(BASE_DIR)/test
OBJ_DIR = $(BASE_DIR)/build
INC_DIR = $(BASE_DIR)/include

PREFIX = $(HOME)

INCLUDES = -I$(HOME)/include -I$(BASE_DIR)/include \
           -I$(BASE_DIR)/lib/include

LFLAGS = -L$(HOME)/lib -L$(BASE_DIR) -L$(BASE_DIR)/lib

LIBS = -lopenblas -lqv
TEST_LIBS = $(LIBS) -lgtest -lgtest_main -lboost_system -pthread

MAIN = $(SRC_DIR)/qvrefl.cc
CARTAN = $(SRC_DIR)/cartan.cc
BENCH = $(SRC_DIR)/benchmark.cc
SRCS = $(filter-out $(MAIN) $(BENCH) $(CARTAN), $(wildcard $(SRC_DIR)/*.cc))
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cc)

_OBJS = $(SRCS:.cc=.o)
_M_OBJ = $(MAIN:.cc=.o)
_B_OBJ = $(BENCH:.cc=.o)
_C_OBJ = $(CARTAN:.cc=.o)
_TEST_OBJS = $(TEST_SRCS:.cc=.o)

OBJS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_OBJS))
M_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_M_OBJ))
B_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_B_OBJ))
C_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_C_OBJ))
TEST_OBJS = $(patsubst $(TEST_DIR)/%,$(OBJ_DIR)/%,$(_TEST_OBJS))

.PHONY: clean all test depend

all:	$(NAME)

$(NAME): $(M_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o $(NAME) $(M_OBJ) $(OBJS) $(LFLAGS) $(LIBS)

$(TEST): $(OBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(B_OPT) $(INCLUDES) -o $(TEST) $(TEST_OBJS) $(OBJS) $(LFLAGS) $(TEST_LIBS)

bench: CXXFLAGS += -Wno-unused-variable
bench: $(B_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o bench $(B_OBJ) $(OBJS) $(LFLAGS) $(LIBS) -lbenchmark -pthread
	@./bench

cartan: $(C_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o cartan $(C_OBJ) $(OBJS) $(LFLAGS) $(LIBS)

test: $(TEST)
	@echo Running tests
	@./$(TEST)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJS): | $(OBJ_DIR)
$(M_OBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	$(RM) *.o *~ $(OBJ_DIR)/*.o $(NAME) $(TEST)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
