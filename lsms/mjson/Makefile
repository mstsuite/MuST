OBJECTS = json.o json_helper.o rstring.o
HEADERS = json.h json_helper.h rstring.h
SRC = json.c json_helper.c rstring.c

all: mjson.a

clean:
	rm -f mjson.a $(OBJECTS)

mjson.a : $(OBJECTS)
	ar -r mjson.a $(OBJECTS)
