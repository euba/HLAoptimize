from setuptools import setup

setup(
    name='HLAopti',
    py_modules=['HLAopti'],
    scripts=['bin/HLAopti_cli','bin/HLAoptiGA_cli','bin/HLAoptiGArand_cli'],
    install_requires=[
        'setuptools',
        'pandas >= 0.22.0',
        'numpy >= 1.16.0'
    ]
)