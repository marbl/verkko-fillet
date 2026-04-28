from setuptools import setup, find_packages

setup(
    name="verkkofillet",                 # Your package name
    version="0.1.21",                      # Version
    author="Juhyun Kim",                   # Author name
    author_email="kimj75@nih.gov", # Author email
    description="A toolkit for cleaning Verkko assemblies.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jjuhyunkim/verkko-fillet", # Repository URL
    packages=find_packages(),             # Automatically find packages in your directory
    # packages=find_packages(include=['verkkofillet', 'verkkofillet.bin', 'verkkofillet.bin.vgp-assembly']),
    # package_dir={'': ''},
    include_package_data=True,
    # package_data={'verkkofillet.bin': ['*sh'], 'verkkofillet.bin.vgp-assembly' : ['*sh']},
    package_data={"":["*"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=[],
)
