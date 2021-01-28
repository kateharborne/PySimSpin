import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PySimSpin", # Replace with your own username
    version="2.0.0",
    author="Katherine Harborne",
    author_email="katherine.harbrone@uwa.edu.au",
    description="A package for producing mock observations of particle simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kateharborne/PySimSpin",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
