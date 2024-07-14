from setuptools import setup

setup(
    name="modena",
    version="1.0.0",
    python_requires=">=3.8",
    package_dir={"": "src"},
    install_requires=[
        "astropy >= 5.2, < 6",
        "click == 8.1.3",
        "h5py == 3.8.0",
        "kmeans1d == 0.3.1",
    ],
    extras_require={'dev': [
        'pytest == 7.2.2',
    ]},
)
