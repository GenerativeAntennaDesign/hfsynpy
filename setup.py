from setuptools import setup, find_packages

setup(
    name="hfsynpy",
    version="0.1.0",
    description="A Python library for synthesis and analysis of high frequency component (currently only microstrip transmission lines), providing accurate models and convenient tools for PCB and RF design.",
    author="Dominik Mair",
    author_email="dominik.mair@uibk.ac.at",
    url="https://github.com/YourUsername/hfsynpy",
    packages=find_packages(),
    install_requires=[
        # No external dependencies required
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
