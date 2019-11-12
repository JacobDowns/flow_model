import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flow_model",
    version="0.1",
    author="Jacob Downs, Doug Brinkerhoff, Jesse Johnson",
    author_email="jacob.downs@umontana.edu",
    description="A flowline ice sheet model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JacobDowns/flow_model",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
