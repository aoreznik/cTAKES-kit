from .definitions import *
import json
from collections import defaultdict, Counter
import pickle
import re
import pyupset as pyu
import pandas as pd
from math import ceil
import hashlib
from warnings import warn
import obonet
from operator import itemgetter


class Element:

    def __repr__(self):
        return "{}: {}".format(str(type(self)), str(self))

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __lt__(self, other):
        return str(self) < str(other)

    def __gt__(self, other):
        return str(self) > str(other)

    def __str__(self):
        raise NotImplementedError


class Disease(Element):

    def __init__(self, id, source, term):
        self.id = id
        self.source = source
        self.term = term

    @property
    def name(self):
        return self.term

    def __str__(self):
        if self.id and self.source:
            return f'{self.source}: {self.id}'
        else:
            return f'No reference: {self.term}'

    _DISEASE_ONTOLOGY_URL = 'https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/HumanDO.obo'
    DISEASE_ONTOLOGY = obonet.read_obo(_DISEASE_ONTOLOGY_URL)


class Drug(Element):

    def __init__(self, id, source, term, **kwargs):
        self.id = id
        self.source = source
        self.term = term

    def __str__(self):
        return str(self.term)


class Gene(Element):

    SYMBOL_TABLE = dict()
    SYMBOL_ALIAS_TABLE = defaultdict(list)
    with open(str(DATA_ROOT / 'non_alt_loci_set.json'), 'r') as f:
        d = json.load(f)
        for doc in d['response']['docs']:
            SYMBOL_TABLE[doc['symbol']] = doc
            for alias in doc.get('alias_symbol', []):
                SYMBOL_ALIAS_TABLE[alias].append(doc['symbol'])
            for prev in doc.get('prev_symbol', []):
                SYMBOL_ALIAS_TABLE[prev].append(doc['symbol'])
        SYMBOL_ALIAS_TABLE = dict(SYMBOL_ALIAS_TABLE)

    def __init__(self, gene_symbol):
        self.gene_symbol = gene_symbol
        try:
            doc = Gene.SYMBOL_TABLE[gene_symbol]
        except KeyError:
            aliases = Gene.SYMBOL_ALIAS_TABLE[gene_symbol]
            # if len(aliases) > 1:
            #     raise KeyError("{} is an ambiguous gene symbol.".format(gene_symbol))
            assert len(aliases) <= 1, 'Ambiguous gene symbol {}'.format(gene_symbol)
            doc = Gene.SYMBOL_TABLE[aliases[0]]
        self.entrez_id = doc['entrez_id']
        self._doc = doc

    @property
    def symbol(self):
        return self.gene_symbol

    def __str__(self):
        return str(self.gene_symbol)

    def __bool__(self):
        return bool(self.entrez_id)

    def __hash__(self):
        return int(self.entrez_id)

    def __eq__(self, other):
        return self.entrez_id == other.entrez_id

    def __lt__(self, other):
        return self.entrez_id < other.entrez_id

    def __gt__(self, other):
        return self.entrez_id > other.entrez_id


class GenomicFeature(Element):

    CHROMOSOMES = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']
    REFERENCE_BUILDS = ['GRCh37', 'GRCh38']

    def __init__(self, chromosome, start, end, referenceName='GRCh37', name='', geneSymbol='', sequence_ontology={}, alt=None, ref=None, **kwargs):
        chromosome = str(chromosome)
        if chromosome.lower().startswith('chr'):
            chromosome = chromosome[3:]
        if chromosome == '23':
            chromosome = 'X'
        if chromosome == '24':
            chromosome = 'Y'
        assert chromosome in GenomicFeature.CHROMOSOMES
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        assert self.start <= self.end
        self.so = sequence_ontology
        self.alt = alt
        self.ref = ref
        self.name = name
        self.gene_symbol = geneSymbol
        assert referenceName in GenomicFeature.REFERENCE_BUILDS
        self.reference_name = referenceName

    def __str__(self):
        return ':'.join([str(getattr(self, x)) for x in ['reference_name', 'chromosome', 'start', 'end', 'name']])

    def __eq__(self, other):
        return all([
            self.chromosome == other.chromosome,
            self.start == other.start,
            self.end == other.end,
            self.reference_name == other.reference_name
        ])

    def __hash__(self):
        return hash(tuple([str(getattr(self, x)) for x in ['reference_name', 'chromosome', 'start', 'end', 'alt']]))

    def issubfeature(self, other):
        return all([
            self.chromosome == other.chromosome,
            self.start >= other.start,
            self.end <= other.end,
            self.reference_name == other.reference_name
        ])

    def issuperfeature(self, other):
        return all([
            self.chromosome == other.chromosome,
            self.start <= other.start,
            self.end >= other.end,
            self.reference_name == other.reference_name
        ])

    def __lt__(self, other):
        if self.reference_name != other.reference_name:
            return self.reference_name < other.reference_name
        elif self.chromosome != other.chromosome:
            c = GenomicFeature.CHROMOSOMES
            return c.index(self.chromosome) < c.index(other.chromosome)
        elif self.start != other.start:
            return self.start < other.start
        elif self.end != other.end:
            return self.end < other.end
        else:
            return False

    def __gt__(self, other):
        return not self < other and self != other

    def __le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __contains__(self, item):
        return self.issuperfeature(item)

    def __len__(self):
        return self.end - self.start + 1


class Publication(Element):

    pmid_re = re.compile(r'https?://.*pubmed/(\d+)$')

    def __init__(self, publication_string):
        pmid_match = Publication.pmid_re.match(publication_string)
        self.pmid = None
        self.publication_string = publication_string
        if pmid_match:
            self.pmid = int(pmid_match[1])

    def __str__(self):
        if self.pmid:
            return str(self.pmid)
        else:
            return self.publication_string


class ViccAssociation(dict):

    def __str__(self):
        return str(hash(self))

    def __hash__(self):
        return self._stable_hash()

    def _stable_hash(self):
        raise NotImplementedError

    @property
    def publications(self):
        evidence = self['association']['evidence']
        all_pubs = list()
        for e in evidence:
            all_pubs += [Publication(p) for p in e['info']['publications'] if p]
        return all_pubs

    @property
    def evidence_level(self):
        return self['association']['evidence_label']

    @property
    def genes(self):
        if getattr(self, '_genes', None):
            return self._genes
        out = list()
        for g in self['genes']:
            if not g:
                continue
            try:
                out.append(Gene(g))
            except KeyError:
                continue
            except AssertionError:
                warn('Ambiguous gene symbol {} in assertion {}'.format(g, self))
                continue
        self._genes = out
        return out

    @property
    def source(self):
        return self['source']

    @property
    def features(self):
        if getattr(self, '_features', None):
            return self._features
        out = list()
        for f in self['features']:
            try:
                f2 = GenomicFeature(**f)
            except:
                continue
            out.append(f2)
        self._features = sorted(out)
        return sorted(out)

    @property
    def disease(self):
        try:
            return Disease(**self['association']['phenotype']['type'])
        except KeyError:
            return None

    @property
    def drugs(self):
        out = list()
        for d in self['association'].get('environmentalContexts', []):
            try:
                d2 = Drug(**d)
            except:
                continue
            out.append(d2)
        return out

    def __eq__(self, other):
        return hash(self) == hash(other)

    @property
    def description(self):
        return self['association'].get('description', None)


class RawAssociation(ViccAssociation):

    def _stable_hash(self):
        source = self['source']
        if source == 'civic':
            assert len(self['association']['evidence']) == 1  # we currently import 1 evidence per association.
            k = 'civic:{}'.format(self['association']['evidence'][0]['evidenceType']['id'])
        elif source == 'molecularmatch':
            k = 'mm:{}'.format(self['raw']['hashKey'])
        elif source == 'brca':
            k = 'brca:{}'.format(self['raw']['id'])
        elif source == 'pmkb':
            t = [self['raw']['variant']['id'], self['raw']['tumor']['id']] + [x['id'] for x in self['raw']['tissues']]
            k = 'pmkb:{}'.format('-'.join([str(x) for x in t]))  # There's no interpretation ID, made compound ID from components
        elif source == 'oncokb':
            try:
                pre = self['raw']['clinical']
                t = [pre['cancerType'], pre['drug'], pre['gene'], pre['variant']['name'], pre['level']]
                k = 'oncokb_clinical:{}'.format('-'.join(t))
            except KeyError:
                pre = self['raw']['biological']
                t = [pre['gene'], pre['variant']['name'], pre['oncogenic'], pre['mutationEffectPmids']]
                k = 'oncokb_biological:{}'.format('-'.join(t))
        elif source == 'jax':
            k = 'jax:{}'.format(self['raw']['id'])
        elif source == 'cgi':
            t = [self['raw']['Drug full name'], self['raw']['Primary Tumor type'],
                 self['raw']['Alteration'], self['raw']['Source']] + self['raw']['individual_mutation']
            k = 'cgi:{}'.format('-'.join(t))
        else:
            raise NotImplementedError("No hash routine defined for source '{}'".format(source))
        b = k.encode()
        m = hashlib.sha256()
        m.update(b)
        return int(m.hexdigest(), 16) % 10**16  # Last 16 digits of sha256 hash are sufficient


class ViccDb:

    DEFAULT_CACHE = DATA_ROOT / 'association_cache.pkl'

    def __init__(self, associations=None, load_cache=False,
                 save_cache=False, cache_path=DEFAULT_CACHE,
                 data_dir=(DATA_ROOT / '0.10')
                 ):
        if load_cache and save_cache:
            raise ValueError('Can only load or save cache, not both.')
        if load_cache:
            with open(str(cache_path), 'rb') as f:
                self.associations = pickle.load(f)
        elif associations is not None:
            self.associations = associations
        else:
            self.load_data(data_dir=data_dir)
        self._index_associations()
        if save_cache:
            self.cache_data(cache_path)

    def load_data(self, data_dir):
        resource_paths = list(data_dir.glob('*.json'))
        if resource_paths:
            self._load_local(resource_paths)
        else:
            self._load_s3()

    def _load_local(self, resource_paths):
        self.associations = list()
        for path in resource_paths:
            source = path.parts[-1].split('.')[0]
            with path.open() as json_data:
                for line in json_data:
                    association = RawAssociation(json.loads(line))  # TODO: Move to ViccAssociation after RawAssociation checks pass
                    association['raw'] = association.pop(source)
                    self.associations.append(association)

    def _load_s3(self):
        raise NotImplementedError

    def _index_associations(self):
        features = []
        associations_by_source = defaultdict(list)
        hashed = defaultdict(list)
        for association in self.associations:
            source = association['source']
            associations_by_source[source].append(association)
            h = hash(association)
            hashed[h].append(association)
            for feature in association.features:
                features.append((feature, h))
        self._features = features
        self._features_sorted = False
        self.associations_by_source = dict(associations_by_source)
        self._hashed = hashed
        self._element_by_source = dict()

    @property
    def features(self):
        if self._features_sorted:
            return self._features
        else:
            self._features = sorted(self._features, key=itemgetter(0))
            self._features_sorted = True
            return self._features

    def select(self, filter_function):
        associations = filter(filter_function, self.associations)
        return ViccDb(list(associations))

    def by_source(self, source):
        return ViccDb(self.associations_by_source[source])

    def report_groups(self, superset=None):
        if superset is None:
            total = len(self)
            for group in sorted(self.associations_by_source):
                count = len(self.associations_by_source[group])
                print("{}: {} ({:.1f}% of total)".format(group, count, count / total * 100))
            print("{} total associations".format(total))
        else:
            for group in sorted(self.associations_by_source):
                count = len(self.associations_by_source[group])
                # intended: below will raise error if key doesn't exit in superset, should be actual superset of self.
                superset_count = len(superset.associations_by_source[group])
                print("{}: {} ({:.1f}% of superset)".format(group, count, count / superset_count * 100))
            print("Total: {} ({:.1f}% of superset)".format(len(self.associations),
                                                           len(self.associations) / len(superset.associations) * 100))

    def cache_data(self, cache_path=DEFAULT_CACHE):
        with open(str(cache_path), 'wb') as f:
            pickle.dump(self.associations, f)

    def __len__(self):
        return len(self.associations)

    def __iter__(self):
        return iter(self.associations)

    def __contains__(self, item):
        return hash(item) in self._hashed

    def __getitem__(self, item):
        return self.associations[item]

    def __sub__(self, other):
        for h, associations in other._hashed.items():
            error_msg = "Cannot perform set substraction, association hash not unique."
            assert len(associations) == 1, error_msg
            assert len(self._hashed.get(h, [])) <= 1, error_msg
            # these assertions assume that hash uniquely identifies an association.
            # Currently not true, but should be with harvester changes.
        return ViccDb([x for x in self if x not in other])

    def plot_element_by_source(self, element, filter_func=lambda x: bool(x), min_bound=1, max_bound=1000000000):
        element_by_source = self.get_element_by_source(element)

        df_dict = dict()
        column_name = ['attribute']
        for source in element_by_source:
            filtered_elements = list(filter(filter_func, element_by_source[source]))
            df_dict[source] = pd.DataFrame(filtered_elements, columns=column_name)
        x = pyu.plot(df_dict, unique_keys=column_name, inters_size_bounds=(min_bound, max_bound))
        x['input_data'] = element_by_source
        return x

    def element_by_source_stats(self, element, filter_func=lambda x: bool(x)):
        element_by_source = self.get_element_by_source(element)
        for source, elements in element_by_source.items():
            element_by_source[source] = set(list(filter(filter_func, elements)))
        ubiquitous_elements = set.intersection(*(element_by_source.values()))
        total_elements = set.union(*(element_by_source.values()))
        count = Counter()
        for source in element_by_source:
            count.update(element_by_source[source])
        majority_size = ceil(len(element_by_source) / 2)
        majority_elements = set([element for element in count if count[element] >= majority_size])
        unique_elements = set([element for element in count if count[element] == 1])
        out = {
            'total': total_elements,
            'ubiquitous': ubiquitous_elements,
            'majority': majority_elements,
            'majority_size': majority_size,
            'unique_elements': unique_elements
        }
        a = len(unique_elements)
        b = len(total_elements)
        print("{} / {} ({:.2%}) of {} are represented in only 1 resource."
              .format(a, b, a / b, element))

        a = len(majority_elements)
        print("{} / {} ({:.2%}) of {} are represented in the majority of ({}) resources."
              .format(a, b, a / b, element, majority_size))

        a = len(ubiquitous_elements)
        print("{} / {} ({:.2%}) of {} are represented across all resources."
              .format(a, b, a/b, element))
        return out

    def get_element_by_source(self, element):
        try:
            e = self._element_by_source[element]
        except KeyError:
            element_by_source = defaultdict(set)
            for association in self:
                association_element = None
                association_element = getattr(association, element)
                if hasattr(association_element, '__iter__') and not isinstance(association_element, str):
                    element_by_source[association.source].update(association_element)
                elif association_element is None:
                    continue
                else:
                    element_by_source[association.source].add(association_element)
            self._element_by_source[element] = dict(element_by_source)
            e = self._element_by_source[element]
        return e

    MATCH_RANKING = ['exact', 'positional', 'focal', 'regional']

    def search_by_feature(self, chromosome=None, start=None, end=None, reference_name=None,
                          name=None, alt=None, gene_symbol=None, genomic_feature=None):
        """
        Returns a list of hits, each corresponding to a single query / association match.
        All features of an association matching a query are stored in the matches attribute.
        The best match between a query and association is stored in the best_match attribute.
        """
        if not isinstance(genomic_feature, GenomicFeature):
            query = GenomicFeature(chromosome, start, end, reference_name, name, gene_symbol, alt=alt)
        else:
            query = genomic_feature
        return self.search_by_features([query])

    @staticmethod
    def _get_best_match(matches):
        if len(matches) == 1:
            best_match = matches[0]
        else:
            s1 = sorted(matches, key=lambda x: x['p'], reverse=True)
            s2 = sorted(s1, key=lambda x: ViccDb.MATCH_RANKING.index(x['type']))
            best_match = s2[0]
        return best_match

    @staticmethod
    def _get_match_type(query, feature):
        match = {'feature': feature}
        if query == feature:
            if query.alt and feature.alt and query.alt == feature.alt:
                match['type'] = 'exact'
            else:
                match['type'] = 'positional'
            match['p'] = 1
        else:
            if query.issubfeature(feature):
                p = len(query) / len(feature)
            elif query.issuperfeature(feature):
                p = len(feature) / len(query)
            else:
                # When query and feature are only partially overlapping
                if query < feature:
                    num = query.end - feature.start + 1
                else:
                    num = feature.end - query.start + 1
                den = max(len(feature), len(query))
                p = num / den
            assert p < 1, f'{p} should be less than 1 for {query} on {feature}'
            match['p'] = p
            if p >= 0.1:
                match['type'] = 'focal'
            elif p > 0:
                match['type'] = 'regional'
            else:
                raise ValueError(f'Expected an overlap between {query} and {feature}')
        return match

    def search_by_features(self, genomic_features):
        assert isinstance(genomic_features, list)
        db_features_pointer = 0
        query_features_pointer = 0
        match_start = None
        last_query_features_pointer = -1
        c = GenomicFeature.CHROMOSOMES
        query_features = sorted(genomic_features)
        hit_index = dict()
        while query_features_pointer < len(query_features) and db_features_pointer < len(self.features):
            if last_query_features_pointer != query_features_pointer:
                q = query_features[query_features_pointer]
                if match_start is not None:
                    db_features_pointer = match_start
                    match_start = None
                last_query_features_pointer = query_features_pointer
            d, association_hash = self.features[db_features_pointer]
            if q.reference_name != d.reference_name:
                raise NotImplementedError('All records in query and datastore currently must match same reference')
            if c.index(q.chromosome) < c.index(d.chromosome):
                query_features_pointer += 1
                continue
            if c.index(q.chromosome) > c.index(d.chromosome):
                db_features_pointer += 1
                continue
            if q.start > d.end:
                db_features_pointer += 1
                continue
            if q.end < d.start:
                query_features_pointer += 1
                continue
            m = ViccDb._get_match_type(q, d)
            key = (q, association_hash)
            matches = hit_index.get(key, list())
            matches.append(m)
            hit_index[key] = matches
            if match_start is None:
                match_start = db_features_pointer
            db_features_pointer += 1
        hits = list()
        for key, matches in hit_index.items():
            q, association_hash = key
            best_match = ViccDb._get_best_match(matches)
            a = self.get_association_by_hash(association_hash)
            hit = {
                'query': q,
                'association': a,
                'matches': matches,
                'best_match': best_match
            }
            hits.append(hit)
        return hits

    def get_association_by_hash(self, key):
        a = self._hashed[key]
        assert len(a) == 1
        return a[0]

    @property
    def sources(self):
        return self.associations_by_source.keys()