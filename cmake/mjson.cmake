add_library(mjson STATIC)

target_sources(mjson PUBLIC
        ${PROJECT_SOURCE_DIR}/external/mjson/json.c
        ${PROJECT_SOURCE_DIR}/external/mjson/json.h
        ${PROJECT_SOURCE_DIR}/external/mjson/json_helper.c
        ${PROJECT_SOURCE_DIR}/external/mjson/json_helper.h
        ${PROJECT_SOURCE_DIR}/external/mjson/rstring.c
        ${PROJECT_SOURCE_DIR}/external/mjson/rstring.h
        )

target_include_directories(
        mjson PUBLIC
        ${PROJECT_SOURCE_DIR}/external/mjson/
)
