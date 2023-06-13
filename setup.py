from setuptools import setup

setup(
    name="nodclust",
    version="0.0.1",
    python_requires=">=3.10",
    package_dir={"": "src"},
    install_requires=[
        "astropy == 5.2.2",
        "click == 8.1.3",
        "h5py == 3.7.0",
        "jenkspy == 0.3.2",
        "scipy == 1.10.0",
    ],
    extras_require={'dev': [
        'pytest == 7.2.2',
    ]},
)
