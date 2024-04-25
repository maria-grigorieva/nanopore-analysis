import re
from fuzzywuzzy import fuzz
from difflib import SequenceMatcher

def find_fuzzy_string_occurrences(input_string, fuzzy_string, threshold=0.8):
    occurrences = []

    for i in range(len(input_string) - len(fuzzy_string) + 1):
        substring = input_string[i:i + len(fuzzy_string)]
        similarity = SequenceMatcher(None, fuzzy_string, substring).ratio()

        if similarity >= threshold:
            occurrences.append(i)

    return occurrences

def barcode_splitter(s, barcode, result):

    indices = []
    left = barcode[:5]
    right = barcode[-5:]
    pattern = rf'{left}(.*?){right}'
    matches = re.finditer(pattern, s['sequence'])
    substrings = [match.group() for match in matches]
    for sub in substrings:
        indices.extend([(m.start(),m.end()) for m in re.finditer(sub, s['sequence'])])

    start = 0
    if len(indices) > 0:
        for i in range(0, len(indices)):
            new_sequence = {'sequence': s['sequence'][start:indices[i][0]],
                           'score': s['score'][start:indices[i][0]]}
            start = indices[i][1]
            if indices[i][0] - start > 51:
                result.append(new_sequence)

        if start < len(s['sequence']):
            if len(s['sequence'])-start > 51:
                result.append({'sequence': s['sequence'][start:],
                               'score': s['score'][start:]})

def barcode_remover(s, barcode, result, mode='both'):
    indices = []
    left = barcode[:6]
    right = barcode[-6:]
    if mode == 'both':
        pattern = rf'{left}(.*?){right}'
        matches = re.finditer(pattern, s['sequence'])
        substrings = [match.group() for match in matches]
        for sub in substrings:
            indices.extend([(m.start(), m.end()) for m in re.finditer(sub, s['sequence'])])

        if len(indices) > 0:
            seq = ""
            score = []
            last_end = 0

            for start, end in indices:
                seq += s['sequence'][last_end:start]
                score.extend(s['score'][last_end:start])
                last_end = end + 1

            seq += s['sequence'][last_end:]
            score.extend(s['score'][last_end:])
            result.append({'sequence': seq, 'score': score})
    elif mode == 'right':
        right = barcode[-8:]
        idx = [m.end() for m in re.finditer(right, s['sequence'])][-1]
        result.append({'sequence': s['sequence'][idx+1:],
                       'score': s['score'][idx+1:]})
    else:
        result.append(s)


def string_splitter(s, left, right, barcode, result):

    split_indices = []

    idx_left = [val for val in [i.start() for i in re.finditer(left, s['sequence'])]]
    idx_right = [val for val in [i.start() for i in re.finditer(right, s['sequence'])]]

    if len(idx_right) > 0 and len(idx_left) > 0:

        for i in idx_left:
            for j in idx_right:
                if j - i >= len(left) and j - i <= len(barcode):
                    split_indices.append((i, j))

        start = 0
        for i in range(0, len(split_indices)):
            if split_indices[i][0] - start > 50:
                result.append({'sequence': s['sequence'][start:split_indices[i][0]],
                               'score': s['score'][start:split_indices[i][0]]})
            start = split_indices[i][1]

        if start < len(s['sequence']):
            if len(s['sequence']) - start+len(left) > 50:
                result.append({'sequence': s['sequence'][start+len(left):],
                               'score': s['score'][start+len(left):]})

    else:
        result.append(s)

def remove_barcodes_from_sequence(s, barcodes, sequences, mode='split'):

    for b in barcodes:
        if mode == 'split':
            barcode_splitter(s, b, sequences)
        else:
            barcode_remover(s, b, sequences, mode='right')
        # string_splitter(s, left, right, b, sequences)
    # remove barcode by the exact match
    # for b in barcodes:
    #     idx = [val for val in [i.start() for i in re.finditer(b, s['sequence'])]]
    #     if len(idx) > 0:
    #         for i in idx:
    #             sequences.append({'sequence': s['sequence'][i + len(b):], 'score': s['score'][i + len(b):]})

    # if len(sequences) == 0:
    #
    #     left = barcodes[0][:8]
    #     right = barcodes[0][-8:]
    #
    #     fuzzy_left = find_fuzzy_string_occurrences(s['sequence'], left)
    #     fuzzy_right = find_fuzzy_string_occurrences(s['sequence'], right)
    #
    #     for i in fuzzy_right:
    #         x = [j for j in fuzzy_left if i - j <= len(barcodes[0])]
    #         if len(x) > 0:
    #             sequences.append(
    #                 {'sequence': s['sequence'][i + len(right):], 'score': s['score'][i + len(right):]})

        #
        # idx_left = [val for val in [i.start() for i in re.finditer(left, s['sequence'])]]
        # idx_right = [val for val in [i.start() for i in re.finditer(right, s['sequence'])]]

        # if len(idx_right) > 0:
        #     idx = max(idx_right)
        #     sequences.append(
        #         {'sequence': s['sequence'][idx+len(right):], 'score': s['score'][idx+len(right):]})

        #
        # for i in range(0, len(idx_left)):
        #     for b in barcodes:
        #         if len(idx_left) > 0:
        #             if fuzz.ratio(s['sequence'][idx_left[i]:idx_left[i] + len(b)], b) >= 0.8:
        #                 sequences.append({'sequence': s['sequence'][:idx_left[i]], 'score': s['score'][:idx_left[i]]})
        #                 sequences.append({'sequence': s['sequence'][idx_left[i] + len(b):], 'score': s['score'][idx_left[i] + len(b):]})
        #
        # if len(sequences) == 0 and len(idx_right) > 0:
        #     for i in range(0, len(idx_right)):
        #         for b in barcodes:
        #             if idx_right[i] - len(b) >= 0:
        #                 if fuzz.ratio(s['sequence'][idx_right[i] + len(right) - len(b):idx_right[i] + len(right)], b) >= 0.8:
        #                     sequences.append(
        #                         {'sequence': s['sequence'][:idx_right[i] + len(right) - len(b)], 'score': s['score'][:idx_right[i] + len(right) - len(b)]})
        #                     sequences.append({'sequence': s['sequence'][idx_right[i] + len(right):],
        #                                       'score': s['score'][idx_right[i] + len(right):]})
        #
        # if len(sequences) == 0:
        #     pattern = rf'{left}(.*?){right}'
        #     matches = re.finditer(pattern, s['sequence'])
        #
        #     substrings = [match.group() for match in matches]
        #     for sub in substrings:
        #         idx = re.finditer(sub, s['sequence'])
        #         sequences.append({'sequence': s['sequence'][:idx], 'score': s['score'][:idx]})
        #         sequences.append({'sequence': s['sequence'][idx+len(sub):], 'score': s['score'][idx+len(sub)]})

    # return sequences
def remove_barcodes_all(seq_list, barcodes):
    modified_sequences = []
    for s in seq_list:
        modified_sequences.append(remove_barcodes_from_sequence(s, barcodes))
    return modified_sequences

def split_sequence_by_barcode(s, barcodes):
    # remove everything before the common right part of barcode
    common_right = barcodes[0][-8:]

    idx = [val for val in [i.start() for i in re.finditer(common_right, s['sequence'])]]
    if len(idx) > 0:
        last_occurrence = idx[-1]
        s['sequence'] = s['sequence'][last_occurrence:]
        s['score'] = s['score'][last_occurrence:]

    return s

def no_barcodes(s, barcodes):
    s = remove_barcodes_from_sequence(s, barcodes)
    s = split_sequence_by_barcode(s, barcodes)
    return s

def sequences_without_barcodes(seq_list, barcodes):
    modified_sequences = []
    for s in seq_list:
        remove_barcodes_from_sequence(s, barcodes, modified_sequences, mode='split')
        # modified_sequences.append(split_sequence_by_barcode(s, barcodes))
    # flattened = [item for sublist in modified_sequences for item in sublist]

    return modified_sequences
