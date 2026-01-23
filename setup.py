"""
Setup script for alarmist package
"""
from setuptools import setup, find_packages

setup(
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
)
