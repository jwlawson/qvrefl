NAME = qvrefl
MAJOR = 0
MINOR = 1
VERSION = $(MAJOR).$(MINOR)

ifeq ($(CXX),g++)
CXXFLAGS += -Wall -Wextra -march=native
#CXXFLAGS += -fno-signed-zeros -fno-math-errno -fno-rounding-math
#CXXFLAGS += -fno-signaling-nans -fno-trapping-math
#CXXFLAGS += -ffinite-math-only
OPT += -O3 -g
else
CXXFLAGS += -Wall -xHOST
OPT += -O3 -ipo
B_OPT += $(OPT)
AR = xiar
endif
CXXFLAGS += -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')

# Using cygwin -std=gnu++11 should be used rather than -std=c++11
ifeq ($(uname_O),Cygwin)
	CXXFLAGS += -std=gnu++11 -DCYGWIN_STOI
endif
ifeq ($(uname_S),Linux)
	CXXFLAGS += -std=c++11 -fPIC
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

LFLAGS = -L$(HOME)/lib -L$(BASE_DIR)/lib

LIBS = -lopenblas -llapack -lqv -lginac
TEST_LIBS = -lgtest -lgtest_main -lopenblas -llapack -pthread \
						-lboost_system

MAIN = $(SRC_DIR)/main.cc
SRCS = $(filter-out $(MAIN), $(wildcard $(SRC_DIR)/*.cc))
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cc)

_OBJS = $(SRCS:.cc=.o)
_M_OBJ = $(MAIN:.cc=.o)
_TEST_OBJS = $(TEST_SRCS:.cc=.o)

OBJS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_OBJS))
M_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_M_OBJ))
TEST_OBJS = $(patsubst $(TEST_DIR)/%,$(OBJ_DIR)/%,$(_TEST_OBJS))

PROF = #-fprofile-use

.PHONY: clean

all:	$(NAME)

$(NAME): $(M_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o $(NAME) $(M_OBJ) $(OBJS) $(LFLAGS) $(LIBS)

$(TEST): $(OBJS) $(TEST_OBJS)
	$(CXX) $(PROF) $(CXXFLAGS) $(B_OPT) $(INCLUDES) -o $(TEST) $(TEST_OBJS) $(OBJS) $(LFLAGS) $(TEST_LIBS)

test: $(TEST)
	@echo Running tests
	@./$(TEST)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(PROF) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(PROF) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cc
	$(CXX) $(PROF) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJS): | $(OBJ_DIR)
$(M_OBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	$(RM) *.o *~ $(OBJ_DIR)/*.o $(NAME) $(TEST)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
