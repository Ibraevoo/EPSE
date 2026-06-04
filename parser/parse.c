#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define DIM 480

typedef struct ProteinEmbedding {
    uint32_t protein_id;
    float embedding[DIM];
} ProteinEmbedding;

typedef struct {
    char *data;
    size_t capacity;
    size_t length;
} sequence;

typedef struct {
    uint32_t id;
    char accession[32];
    char entry_name[32];
    char description[256];
    char organism[256];
    sequence *seq;
} Entry;

// Utilities for sequence dynamic string container
void sequence_init(sequence *s) {
    s->capacity = 1024;
    s->length = 0;
    s->data = malloc(s->capacity);
    if (s->data == NULL) {
        printf("Malloc failed\n");
        exit(1);
    }
    s->data[0] = '\0';
}

void sequence_append(sequence *s, const char *text) {
    size_t text_len = strlen(text);
    while (s->length + text_len + 1 > s->capacity) {
        s->capacity *= 2;
        char *tmp = realloc(s->data, s->capacity);
        if (tmp == NULL) {
            printf("Realloc failed\n");
            exit(1);
        }
        s->data = tmp;
    }
    memcpy(s->data + s->length, text, text_len);
    s->length += text_len;
    s->data[s->length] = '\0';
}

void sequence_clear(sequence *s) {
    s->length = 0;
    s->data[0] = '\0';
}

void sequence_free(sequence *s) {
    free(s->data);

    s->data = NULL;
    s->length = 0;
    s->capacity = 0;
}

// main function starts here -
void parse_and_save(char *filepath) {

    if (filepath == NULL) {
        printf("\033[31;1;4mERROR: NULL filepath \n");
        return;
    }

    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("\033[31;1;4mERROR: %s could not be opened \n", filepath);
        return;
    }

    FILE *csv = fopen(".data/proteins.csv", "w");

    if (csv == NULL) {
        printf("Failed to create proteins.csv\n");
        free(fp);
        return;
    }

    fprintf(csv, "id,accession,entry_name,description,organism,seq_len\n");

    // parsing logic

    char BUFFER[4096];

    Entry *entry = malloc(sizeof(Entry));
    ProteinEmbedding *embedding = malloc(sizeof(ProteinEmbedding));
    if (entry == NULL || embedding == NULL) {
        printf("m\033[31;1;4mMalloc failed\n");
        fclose(csv);
        fclose(fp);
        return;
    }

    uint32_t id = 0;
    sequence *s = malloc(sizeof(sequence));
    if (s == NULL) {
        printf("m\033[31;1;4mMalloc failed\n");
        fclose(csv);
        fclose(fp);
        free(entry);
        free(embedding);
        return;
    }

    int first_protein = 1;

    sequence_init(s);

    while (fgets(BUFFER, sizeof(BUFFER), fp) != NULL) { // Buffer has the line
        if (BUFFER[0] == '>') {

            if (!first_protein) {

                entry->seq = s;

                fprintf(csv, "%u,\"%s\",\"%s\",\"%s\",\"%s\",%zu\n", entry->id, entry->accession,
                        entry->entry_name, entry->description, entry->organism, s->length);

                printf("Saved protein %u (%s) length=%zu\n", entry->id, entry->accession,
                       s->length);

                sequence_clear(s);
            }

            first_protein = 0;
            id++;
            entry->id = id;

            // New protein header - parse from start to end (end is denoted if we get
            // a space char) - skip 3 char sp| then u get ur id which is till the next
            // |
            int ptr = 4;
            char accession[32];
            int i = 0;
            while (BUFFER[ptr] != '|') {
                accession[i] = BUFFER[ptr];
                i++;
                ptr++;
            }
            accession[i] = '\0';
            // Then parse ascession  skip one char in ptr as its currently '|'

            ptr++; // now it begins at head of the ascession

            char entry_name[32];
            i = 0;
            while (BUFFER[ptr] != ' ') {
                entry_name[i] = BUFFER[ptr];
                i++;
                ptr++;
            }
            entry_name[i] = '\0';
            // then name skip one char in ptr as it points to ' '

            ptr++; // now points to head of description

            char descp[256];

            i = 0;
            while (!(BUFFER[ptr] == 'O' && BUFFER[ptr + 1] == 'S' &&
                     BUFFER[ptr + 2] == '=')) { // get to O of OR

                descp[i] = BUFFER[ptr];
                ptr++;
                i++;
            }
            descp[i] = '\0';
            // ptr now points to O of OR make it point to the first letter of organism
            // name via skipping R and =
            ptr += 3; // goto R goto = now at beginning
            char org[256];
            i = 0;
            while (!(BUFFER[ptr] == 'O' && BUFFER[ptr + 1] == 'X' && BUFFER[ptr + 2] == '=')) {

                org[i] = BUFFER[ptr];
                ptr++;
                i++;
            }
            org[i] = '\0';

            strcpy(entry->accession, accession);
            strcpy(entry->entry_name, entry_name);
            strcpy(entry->description, descp);
            strcpy(entry->organism, org);

        } else {
            // Sequence line

            // Append amino acids to current sequence
            BUFFER[strcspn(BUFFER, "\n")] = '\0';
            sequence_append(s, BUFFER);
        }
    }

    if (!first_protein) {
        fprintf(csv, "%u,\"%s\",\"%s\",\"%s\",\"%s\",%zu\n", entry->id, entry->accession,
                entry->entry_name, entry->description, entry->organism, s->length);
    }

    sequence_free(s);

    free(s);
    free(entry);
    free(embedding);

    fclose(csv);
    fclose(fp);
}

int main() {
    mkdir(".data", 0755);
    parse_and_save("uniprot_sprot.fasta");
}
