import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='decon2vcf',
    version='1.0.0',
    author='AWMGS',
    author_email='bioinformatics.team@nhs.wales.uk',
    description='Convert Decon Output to VCF format.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/AWGL/DeconToVCF',
    packages=setuptools.find_packages(),
    scripts= ['Decon2VCF.py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    install_requires=[
   'pandas>=0.23.4'
],
)