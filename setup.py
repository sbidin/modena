from setuptools import setup

setup(
    name="signals",
    version="0.0.1",
    python_requires=">=3.10",
    package_dir={"": "src"},
    install_requires=[
        "click == 8.*",
        "astropy == 5.*",
        "h5py == 3.*",
    ],
)
