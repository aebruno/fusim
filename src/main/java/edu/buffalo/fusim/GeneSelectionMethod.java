/*
 * Copyright 2012 Andrew E. Bruno <aebruno2@buffalo.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

package edu.buffalo.fusim;

/**
 * Enum for encapsulating Gene selection methods
 * 
 * @author Andrew E. Bruno
 * 
 */
public enum GeneSelectionMethod {
    UNIFORM {
        public String toString() {
            return "uniform";
        }
    },
    BINNED {
        public String toString() {
            return "binned";
        }
    },
    EMPIRICAL_STURGES {
        public String toString() {
            return "empirical-sturges";
        }
    },
    EMPIRICAL {
        public String toString() {
            return "empirical";
        }
    };
    
    public static GeneSelectionMethod fromString(String str) {
        if(str.equalsIgnoreCase(GeneSelectionMethod.UNIFORM.toString())) {
            return GeneSelectionMethod.UNIFORM;
        } else if(str.equalsIgnoreCase(GeneSelectionMethod.BINNED.toString())) {
            return GeneSelectionMethod.BINNED;
        } else if(str.equalsIgnoreCase(GeneSelectionMethod.EMPIRICAL.toString())) {
            return GeneSelectionMethod.EMPIRICAL;
        } else if(str.equalsIgnoreCase(GeneSelectionMethod.EMPIRICAL_STURGES.toString())) {
            return GeneSelectionMethod.EMPIRICAL_STURGES;
        } else {
            return null;
        }
    }
}
