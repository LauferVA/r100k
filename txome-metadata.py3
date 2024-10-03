from Bio import Entrez
import time
import csv
import re

# Set up Entrez email and api_key
Entrez.email = "your_email@example.com"
Entrez.api_key = "YOUR_API_KEY"

# Function to chunk a list into batches
def batcher(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

# Search GEO for whole transcriptome studies
def search_geo_whole_transcriptome():
    search_term = '(whole transcriptome[All Fields]) AND "Expression profiling by high throughput sequencing"[DataSet Type]'
    # First, perform the search to get the total number of results
    handle = Entrez.esearch(db="gds", term=search_term, retmax=0)
    search_results = Entrez.read(handle)
    handle.close()
    count = int(search_results['Count'])
    print(f"Found {count} records.")
    id_list = []
    batch_size = 10000
    for start in range(0, count, batch_size):
        print(f"Fetching IDs {start} to {start+batch_size}")
        handle = Entrez.esearch(db="gds", term=search_term, retstart=start, retmax=batch_size)
        search_results = Entrez.read(handle)
        handle.close()
        id_list.extend(search_results['IdList'])
        time.sleep(0.1)  # Sleep to avoid exceeding rate limits
    return id_list

# Extract accession numbers from text
def extract_accessions(text, prefix):
    pattern = prefix + r'\d+'
    return re.findall(pattern, text)

# Fetch metadata from BioSample
def fetch_biosample_metadata(accession):
    # Query the BioSample database for metadata
    try:
        handle = Entrez.esummary(db='biosample', id=accession)
        summary = Entrez.read(handle)
        handle.close()
        # Extract relevant fields
        metadata = {}
        if 'DocumentSummarySet' in summary:
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            metadata['BioSampleAccession'] = docsum.get('Accession', '')
            metadata['BioSampleTitle'] = docsum.get('Title', '')
            # Add more fields as needed
        return metadata
    except Exception as e:
        print(f"Error fetching BioSample metadata for {accession}: {e}")
        return {}

# Fetch metadata from BioProject
def fetch_bioproject_metadata(accession):
    # Query the BioProject database for metadata
    try:
        handle = Entrez.esummary(db='bioproject', id=accession)
        summary = Entrez.read(handle)
        handle.close()
        # Extract relevant fields
        metadata = {}
        if 'DocumentSummarySet' in summary:
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            metadata['BioProjectAccession'] = docsum.get('Project_Acc', '')
            metadata['BioProjectTitle'] = docsum.get('Project_Title', '')
            # Add more fields as needed
        return metadata
    except Exception as e:
        print(f"Error fetching BioProject metadata for {accession}: {e}")
        return {}

# Extract metadata from a summary record
def extract_metadata(summary):
    record = {}
    # Extract fields from the summary record
    record['ID'] = summary['Id']
    record['Title'] = summary.get('title', '')
    record['Organism'] = summary.get('taxon', '')
    record['Platform'] = summary.get('GPL', '')
    record['Samples'] = summary.get('n_samples', '')
    record['Summary'] = summary.get('summary', '')
    record['Accession'] = summary.get('Accession', '')
    record['PubMedIDs'] = ','.join(summary.get('PubMedIds', []))
    record['SupplementaryData'] = ','.join(summary.get('supplementary_data', []))
    record['SubmissionDate'] = summary.get('SubmissionDate', '')
    record['UpdateDate'] = summary.get('UpdateDate', '')
    record['DataSets'] = ','.join(summary.get('GDS', []))
    record['Series'] = ','.join(summary.get('GSE', []))
    record['SamplesList'] = ','.join(summary.get('GSM', []))
    record['StudyType'] = summary.get('gdstype', '')
    return record

# Cross-reference with BioSample and BioProject
def cross_reference_single_record(record):
    # Extract BioProject or BioSample accession numbers from record
    biosample_accessions = extract_accessions(record.get('Summary', ''), prefix='SAMN')
    bioproject_accessions = extract_accessions(record.get('Summary', ''), prefix='PRJNA')
    # Query BioSample
    for accession in biosample_accessions:
        biosample_metadata = fetch_biosample_metadata(accession)
        record.update(biosample_metadata)
    # Query BioProject
    for accession in bioproject_accessions:
        bioproject_metadata = fetch_bioproject_metadata(accession)
        record.update(bioproject_metadata)
    return record

# Fetch metadata and write to file
def fetch_and_write_metadata(id_list, filename):
    fieldnames = None
    with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = None
        for id_batch in batcher(id_list, 5000):
            ids = ','.join(id_batch)
            print(f"Fetching metadata for IDs: {id_batch[0]} to {id_batch[-1]}")
            handle = Entrez.esummary(db="gds", id=ids)
            summaries = Entrez.read(handle)
            handle.close()
            records = []
            for summary in summaries['DocumentSummarySet']['DocumentSummary']:
                record = extract_metadata(summary)
                # Cross-reference
                record = cross_reference_single_record(record)
                records.append(record)
            # Write records to file
            if writer is None:
                # Get all fieldnames
                fieldnames = set()
                for record in records:
                    fieldnames.update(record.keys())
                fieldnames = list(fieldnames)
                writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
                writer.writeheader()
            for record in records:
                writer.writerow(record)
            time.sleep(0.1)  # Sleep to avoid exceeding rate limits

def main():
    id_list = search_geo_whole_transcriptome()
    fetch_and_write_metadata(id_list, 'geo_whole_transcriptome_metadata.tsv')

if __name__ == '__main__':
    main()
