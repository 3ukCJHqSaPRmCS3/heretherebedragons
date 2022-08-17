import setuptools

with open("README.md", "r") as rfh:
    long_project_desc = rfh.read()

setuptools.setup(
    name="swave", # Sphereflow wave
    version="0.1.1",
    packages=setuptools.find_packages(include=["swave","swave*"]),
    long_description=long_project_desc,
    long_description_content_type="text/markdown",
    install_requires=[
        "cfs>=0.0.2",
        "matplotlib>=3.4",
        "numba>=0.52",
        # "numba>=0.54.0",
        "numpy>=1.19",
        "pandas>=1.3.4",
        "path",
        "pyrsistent>=0.18",
        # "sympy>=1.9"        
    ],
    package_data={
        "propagator":[
            "base/data/*.pkl", #any and all pickle objects deposited in the "data" folder
            "base/data/*.csv" #any and all csv's in the "data" folder
        ]
    }
)
