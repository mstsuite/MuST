all: wltool

wltool: wltool.cpp
	$(CXX) -o wltool -O3 wltool.cpp ../../libjson/libjson.a
	cp wltool $(TOP_DIR)/bin/wltool

clean:
	rm -f wltool $(TOP_DIR)/bin/wltool
