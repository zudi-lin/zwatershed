def setup_cython():
	from distutils.core import setup
	from distutils.extension import Extension
	from Cython.Distutils import build_ext

	ext_modules = [Extension("zwatershed",
                         ["zwatershed.pyx", "zwatershed_main.cpp"],
                         library_dirs=["."],
                         include_dirs=[".",'/n/home04/donglai/.conda/envs/zwatershed/include/'],
                         language='c++',
                         # std= 'c++11',
                         extra_link_args=["-std=c++11"],
                         extra_compile_args=["-std=c++11", "-w"])]

	setup(cmdclass={'build_ext': build_ext}, ext_modules=ext_modules)
if __name__=='__main__':
	setup_cython()
