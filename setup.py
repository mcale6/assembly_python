from setuptools import setup, find_packages

setup(
    name="assembly_python",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "rdkit>=2023.3",
        "numpy",
    ],
    entry_points={
        "console_scripts": [
            "assembly=assembly_python.assembly:main",
        ],
    },
    author="mcale6",
    description="Assembly Index calculation in Python",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
)