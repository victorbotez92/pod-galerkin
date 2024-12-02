set(CMAKE_CXX_COMPILER "/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/bin/intel64/icpc")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "19.1.3.20200925")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "98")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "GNU")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_CXX_SIMULATE_VERSION "4.8.5")



set(CMAKE_AR "/bin/ar")
set(CMAKE_CXX_COMPILER_AR "")
set(CMAKE_RANLIB "/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "")
set(CMAKE_LINKER "/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/pstl/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/python-3.9.10-lt3tzy6xiqysondaoxripdbydpruw2wi/include/python3.9;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-mpi-2019.9.304-d4pbusne53goxmkp4drzqim2rbogfmit/compilers_and_libraries/linux/mpi/intel64/include;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/parmetis-4.0.3-hh7ms6hxlc62pcukm37utslljfovv7yi/include;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/intel-mpi-2019.3.199-rplhb3fyzslp2dk6m5i3cwcoc73roagj/compilers_and_libraries/linux/mpi/intel64/include;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/metis-5.1.0-vzvpikojvj3tcczmdfqql6enuxncwief/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/ipp/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/mkl/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/pstl/stdlib;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/tbb/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/daal/include;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/compiler/include/intel64;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/compiler/include/icc;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/compiler/include;/usr/include/c++/4.8.5;/usr/include/c++/4.8.5/x86_64-redhat-linux;/usr/include/c++/4.8.5/backward;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "imf;svml;irng;stdc++;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-mpi-2019.9.304-d4pbusne53goxmkp4drzqim2rbogfmit/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib;/gpfs/softs/libraries/FFTW/3.3.10/lib;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/parmetis-4.0.3-hh7ms6hxlc62pcukm37utslljfovv7yi/lib;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/intel-mpi-2019.3.199-rplhb3fyzslp2dk6m5i3cwcoc73roagj/compilers_and_libraries_2019.3.199/linux/mpi/intel64/libfabric/lib;/gpfs/softs/spack/opt/spack/linux-centos7-cascadelake/intel-19.0.3.199/metis-5.1.0-vzvpikojvj3tcczmdfqql6enuxncwief/lib;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.4;/gpfs/softs/spack_0.17/opt/spack/linux-centos7-cascadelake/intel-20.0.4.304/intel-parallel-studio-cluster.2020.4-o46drv5pkuebzwnoavrfodb5x7ggnzeg/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.8;/usr/lib/gcc/x86_64-redhat-linux/4.8.5;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
