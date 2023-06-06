from setuptools import setup

setup(
    name="nodclust",
    version="0.0.1",
    python_requires=">=3.10",
    package_dir={"": "src"},
    install_requires=[
        "astropy == 5.*",
        "click == 8.*",
        "h5py == 3.*",
        "jenkspy == 0.*",
        "scipy == 1.*",
    ],
    extras_require={'dev': [
        'pytest == 7.*',
    ]},
)
