package com.pescamillam.fem.model;

/**
 * Pojo class with the input values that the program can receive
 * 
 * @author Peter Escamilla (pescamilla@unab.edu.co)
 */
public class InputValues {
    private String width;
    private String height;
    private String thickness;
    private String area;
    private String elasticity;
    private String density;
    private String poisson;
    private String deltaTime;
    private String numTimes;
    private String numX;
    private String numY;
    private String initialDisplacement;

    public InputValues(String width, String height, String thickness, String area, String elasticity, String density,
            String poisson, String deltaTime, String numTimes, String numX, String numY, String initialDisplacement) {
        this.width = width;
        this.height = height;
        this.thickness = thickness;
        this.area = area;
        this.elasticity = elasticity;
        this.density = density;
        this.poisson = poisson;
        this.deltaTime = deltaTime;
        this.numTimes = numTimes;
        this.numX = numX;
        this.numY = numY;
        this.initialDisplacement = initialDisplacement;
    }

    public String getWidth() {
        return width;
    }

    public String getHeight() {
        return height;
    }

    public String getThickness() {
        return thickness;
    }



    public String getArea() {
        return area;
    }



    public String getElasticity() {
        return elasticity;
    }



    public String getDensity() {
        return density;
    }



    public String getPoisson() {
        return poisson;
    }



    public String getDeltaTime() {
        return deltaTime;
    }



    public String getNumTimes() {
        return numTimes;
    }

    public String getNumX() {
        return numX;
    }

    public String getNumY() {
        return numY;
    }

    public String getInitialDisplacement() {
        return initialDisplacement;
    }
    /**
     * Builder to create a InputValues instance easily 
     * 
     * @author Peter Escamilla (pescamilla@unab.edu.co)
     */
    public static class Builder {
        private String width;
        private String height;
        private String thickness;
        private String area;
        private String elasticity;
        private String density;
        private String poisson;
        private String deltaTime;
        private String numTimes;
        private String numX;
        private String numY;
        private String initialDisplacement;

        public Builder withWidth(String width) {
            this.width = width;
            return this;
        }

        public Builder withHeight(String height) {
            this.height = height;
            return this;
        }
        
        public Builder withThickness(String thickness) {
            this.thickness = thickness;
            return this;
        }
        
        public Builder withArea(String area) {
            this.area = area;
            return this;
        }
        
        public Builder withElasticity(String elasticity) {
            this.elasticity = elasticity;
            return this;
        }
        
        public Builder withDensity(String density) {
            this.density = density;
            return this;
        }
        
        public Builder withPoisson(String poisson) {
            this.poisson = poisson;
            return this;
        }
        
        public Builder withDeltaTime(String deltaTime) {
            this.deltaTime = deltaTime;
            return this;
        }
        
        public Builder withNumTimes(String numTimes) {
            this.numTimes = numTimes;
            return this;
        }
        
        public Builder withNumX(String numX) {
            this.numX = numX;
            return this;
        }
        
        public Builder withNumY(String numY) {
            this.numY = numY;
            return this;
        }

        public Builder withInitialDisplacement(String initialDisplacement) {
            this.initialDisplacement = initialDisplacement;
            return this;
        }
        
        public InputValues build() {
            return new InputValues(width, height, thickness, area, elasticity, density, poisson,
                    deltaTime, numTimes, numX, numY, initialDisplacement);
        }
        
    }
}
