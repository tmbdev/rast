// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

/// \file classifier.h
/// \brief Interfaces for classification and density estimation

#ifndef h_classifier__
#define h_classifier__

#include <stdio.h>
#include "colib/narray.h"
#include "colib/narray-util.h"
#include "colib/rowarrays.h"
#include "smartptr.h"

namespace colib {
    struct IClassifier {
        virtual ~IClassifier() {}

        IClassifier() {
            verbose = 0;
        }

        // Number of input features expected.

        virtual int nfeatures() { 
            throw "unimplemented"; 
        }

        // Number of possible output classes.

        virtual int nclasses() { 
            throw "unimplemented"; 
        }

        // Complexity parameter of the classifier (if any).

        virtual float complexity() { 
            return 1.0; 
        }

        // Parameter settings.

        int verbose;

        virtual void set(const char *name,double value) {
            if(!strcmp("verbose",name)) verbose = (int)value;
            else throw "unknown parameter";
        }

        // Various interfaces to classification; there are default
        // implementations in terms of posteriors, discriminants, etc.
        // Override whatever functionality you can actually provide.

        virtual int classify(floatarray &v) {
            floatarray temp;
            discriminant(temp,v);
            return argmax(temp);
        }

        virtual void discriminant(floatarray &result,floatarray &v)  {
            posterior(result,v);
        }

        virtual void posterior(floatarray &result,floatarray &v)  { 
            throw "unimplemented"; 
        }

        virtual void transform(floatarray &out,floatarray &v) {
            throw "unimplemented"; 
        }

        // Negative log likelihood of the sample.

        virtual double cost(floatarray &v) { 
            throw "unimplemented"; 
        }

        // Combination methods computing negative log likelihood and posterior
        // simultaneously (they are separate estimates and not necessarily
        // consistent).

        virtual int classify(float &cost,floatarray &v) {
            cost = this->cost(v);
            return classify(v);
        }

        virtual void discriminant(float &cost,floatarray &result,floatarray &v) {
            cost = this->cost(v);
            discriminant(result,v);
        }

        virtual void posterior(float &cost,floatarray &result,floatarray &v) { 
            cost = this->cost(v);
            posterior(result,v);
        }

        // Batch training.

        virtual void train(floatarray &data,intarray &classes) { 
            throw "unimplemented"; 
        }

        virtual void train(floatarray &data) { 
            throw "unimplemented"; 
        }

        virtual void startTraining() {
               throw "unimplemented";
        }
        virtual void finishTraining() {
               throw "unimplemented";
        }

        // Incremental training and its default implementation.
        // Override these methods if you have better incremental
        // training available for your particular class.
        // Both supervised and unsupervised training methods are
        // provided.  Override and/or throw errors as necessary.

        floatarray vectors;
        intarray classes;

        virtual void add(floatarray &v,int c) {
            rowpush(vectors,v);
            classes.push(c);
        }

        floatarray unsupervised;

        virtual void add(floatarray &v) {
            rowpush(unsupervised,v);
        }

        virtual void train() {
            train(unsupervised);
            train(vectors,classes);
        }

        // Saving and loading the classifier.

        virtual void save(FILE *stream) { 
            throw "unimplemented"; 
        }

        virtual void load(FILE *stream) { 
            throw "unimplemented"; 
        }

        virtual void save(const char *path) {
            save(stdio(path,"wb"));
        }

        virtual void load(const char *path) {
            load(stdio(path,"rb"));
        }

        virtual void print() {
            printf("<some classifier %p>\n",this);
        }
    };
}

#endif
