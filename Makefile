NAME = qvrefl
MAJOR = 0
MINOR = 1
VERSION = $(MAJOR).$(MINOR)

ifeq ($(CXX),icpc)
CXXFLAGS = -Wall -xHOST
OPT = -O3 -ipo -no-prec-div
AR = xiar
endif
ifeq ($(CXX),clang++)
CXXFLAGS += -Wall -Wextra -Werror -march=native -Wno-unused-parameter\
CXXFLAGS += -fno-signed-zeros -fno-math-errno
CXXFLAGS += -fno-trapping-math
CXXFLAGS += -ffinite-math-only
OPT = -Ofast
AR = ar
endif
ifeq ($(CXX),g++)
CXXFLAGS += -Wall -Wextra -march=native
CXXFLAGS += -fno-signed-zeros -fno-math-errno -fno-rounding-math
CXXFLAGS += -fno-signaling-nans -fno-trapping-math
CXXFLAGS += -ffinite-math-only -Wno-misleading-indentation
OPT = -g -Ofast
AR = gcc-ar
endif
override CXXFLAGS += -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -DNDEBUG

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')

ifeq ($(uname_O),Cygwin)
override CXXFLAGS += -std=gnu++1z -DCYGWIN_STOI
endif
ifeq ($(uname_S),Linux)
override CXXFLAGS += -std=c++1z
endif

TEST = test$(NAME)

LIB = lib$(NAME).so.$(VERSION)
STATIC = lib$(NAME).a
LDFLAGS = -shared -Wl,-soname,$(LIB)

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

BINS = qvrefl cartan cmut cexch
MAIN = $(SRC_DIR)/qvrefl.cc
CARTAN = $(SRC_DIR)/cartan.cc
BENCH = $(SRC_DIR)/benchmark.cc
CMUT = $(SRC_DIR)/cmut.cc
CEXCH = $(SRC_DIR)/cexch.cc
SRCS = $(filter-out $(MAIN) $(BENCH) $(CARTAN) $(CMUT) $(CEXCH), $(wildcard $(SRC_DIR)/*.cc))
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cc)

_OBJS = $(SRCS:.cc=.o)
_M_OBJ = $(MAIN:.cc=.o)
_B_OBJ = $(BENCH:.cc=.o)
_C_OBJ = $(CARTAN:.cc=.o)
_CM_OBJ = $(CMUT:.cc=.o)
_CE_OBJ = $(CEXCH:.cc=.o)
_TEST_OBJS = $(TEST_SRCS:.cc=.o)

OBJS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_OBJS))
M_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_M_OBJ))
B_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_B_OBJ))
C_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_C_OBJ))
CM_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_CM_OBJ))
CE_OBJ = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_CE_OBJ))
TEST_OBJS = $(patsubst $(TEST_DIR)/%,$(OBJ_DIR)/%,$(_TEST_OBJS))

.PHONY: clean all test depend

all:	$(BINS)

lib:	$(LIB)
static:	$(STATIC)

$(LIB): CXXFLAGS+=-fPIC
$(LIB):	$(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(STATIC): OPT += -flto -fuse-linker-plugin -ffat-lto-objects
$(STATIC): $(OBJS)
	$(AR) rcs $(STATIC) $(OBJS)

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

cmut: $(CM_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o cmut $(CM_OBJ) $(OBJS) $(LFLAGS) $(LIBS)

cexch: $(CE_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDES) -o cexch $(CE_OBJ) $(OBJS) $(LFLAGS) $(LIBS)

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

install:	$(LIB)
	cp $(LIB) $(PREFIX)/lib
	ldconfig -v -n $(PREFIX)/lib
	ln -fs $(PREFIX)/lib/$(LIB) $(PREFIX)/lib/lib$(NAME).so
	mkdir -p $(PREFIX)/include/$(NAME)
	cp -ru include/* $(PREFIX)/include/$(NAME)/

uninstall:
	rm $(PREFIX)/lib/lib$(NAME).*
	rm -r $(PREFIX)/include/$(NAME)

clean:
	$(RM) *.o *~ $(OBJ_DIR)/*.o $(TEST) bench $(BINS)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
