include(FindPackageHandleStandardArgs)

# search paths
set(GLM_SEARCH_PATH /usr/include /usr/local/include)

# locate include path
find_path(GLM_INCLUDE_DIR "glm/glm.hpp" PATHS ${GLM_SEARCH_PATH})
find_package_handle_standard_args(GLM DEFAULT_MSG GLM_INCLUDE_DIR)

if(GLM_FOUND)
set(GLM_INCLUDE_DIRS "${GLM_INCLUDE_DIR}")
endif(GLM_FOUND)
