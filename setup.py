import setuptools

setuptools.setup(
    name = "gdc_rnaseq_tools",
    author = "Kyle Hernandez",
    author_email = "kmhernan@uchicago.edu",
    version = 0.5,
    description = "helper utilities for the GDC RNA-Seq workflow",
    url = "https://github.com/NCI-GDC/gdc-rnaseq-tool",
    license = "Apache 2.0",
    packages = setuptools.find_packages(),
    entry_points = {
        'console_scripts': [
            'gdc-rnaseq-tools = gdc_rnaseq_tools.__main__:main'
        ]
    }
)
