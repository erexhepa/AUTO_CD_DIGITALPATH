//package ij.plugin.filter;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Random;
////package com.dataonfocus.clustering;
//
// public class Cluster {
//
//        public List points;
//        public Point centroid;
//        public int id;
//
//        //Creates a new Cluster
//        public Cluster(int id) {
//            this.id = id;
//            this.points = new ArrayList();
//            this.centroid = null;
//        }
//
//        public List getPoints() {
//            return points;
//        }
//
//        public void addPoint(Point point) {
//            points.add(point);
//        }
//
//        public void setPoints(List points) {
//            this.points = points;
//        }
//
//        public Point getCentroid() {
//            return centroid;
//        }
//
//        public void setCentroid(Point centroid) {
//            this.centroid = centroid;
//        }
//
//        public int getId() {
//            return id;
//        }
//
//        public void clear() {
//            points.clear();
//        }
//
//        public void plotCluster() {
//            System.out.println("[Cluster: " + id+"]");
//            System.out.println("[Centroid: " + centroid + "]");
//            System.out.println("[Points: \n");
//            for(Point p : points) {
//                System.out.println(p);
//            }
//            System.out.println("]");
//        }
//
//    }
//
//
//public class Point {
//
//
//    private double x = 0;
//    private double y = 0;
//    private int cluster_number = 0;
//
//    public Point(double x, double y)
//    {
//        this.setX(x);
//        this.setY(y);
//    }
//
//    public void setX(double x) {
//        this.x = x;
//    }
//
//    public double getX()  {
//        return this.x;
//    }
//
//    public void setY(double y) {
//        this.y = y;
//    }
//
//    public double getY() {
//        return this.y;
//    }
//
//    public void setCluster(int n) {
//        this.cluster_number = n;
//    }
//
//    public int getCluster() {
//        return this.cluster_number;
//    }
//
//    //Calculates the distance between two points.
//    protected static double distance(Point p, Point centroid) {
//        return Math.sqrt(Math.pow((centroid.getY() - p.getY()), 2) + Math.pow((centroid.getX() - p.getX()), 2));
//    }
//
//    //Creates random point
//    protected static Point createRandomPoint(int min, int max) {
//        Random r = new Random();
//        double x = min + (max - min) * r.nextDouble();
//        double y = min + (max - min) * r.nextDouble();
//        return new Point(x,y);
//    }
//
//    protected static List createRandomPoints(int min, int max, int number) {
//        List points = new ArrayList(number);
//        for(int i = 0; i &lt; number; i++) {
//            points.add(createRandomPoint(min,max));
//        }
//        return points;
//    }
//
//    public String toString() {
//        return "("+x+","+y+")";
//    }
//}
//
//public class KMeans {
//
//    //Number of Clusters. This metric should be related to the number of points
//    private int NUM_CLUSTERS = 3;
//    //Number of Points
//    private int NUM_POINTS = 15;
//    //Min and Max X and Y
//    private static final int MIN_COORDINATE = 0;
//    private static final int MAX_COORDINATE = 10;
//
//    private List points;
//    private List clusters;
//
//    public KMeans() {
//        this.points = new ArrayList();
//        this.clusters = new ArrayList();
//    }
//
//    public static void main(String[] args) {
//
//        KMeans kmeans = new KMeans();
//        kmeans.init();
//        kmeans.calculate();
//    }
//
//    //Initializes the process
//    public void init() {
//        //Create Points
//        points = Point.createRandomPoints(MIN_COORDINATE,MAX_COORDINATE,NUM_POINTS);
//
//        //Create Clusters
//        //Set Random Centroids
//        for (int i = 0; i &lt; NUM_CLUSTERS; i++) {
//            Cluster cluster = new Cluster(i);
//            Point centroid = Point.createRandomPoint(MIN_COORDINATE,MAX_COORDINATE);
//            cluster.setCentroid(centroid);
//            clusters.add(cluster);
//        }
//
//        //Print Initial state
//        plotClusters();
//    }
//
//    private void plotClusters() {
//        for (int i = 0; i &lt; NUM_CLUSTERS; i++) {
//            Cluster c = clusters.get(i);
//            c.plotCluster();
//        }
//    }
//
//    //The process to calculate the K Means, with iterating method.
//    public void calculate() {
//        boolean finish = false;
//        int iteration = 0;
//
//        // Add in new data, one at a time, recalculating centroids with each new one.
//        while(!finish) {
//            //Clear cluster state
//            clearClusters();
//
//            List lastCentroids = getCentroids();
//
//            //Assign points to the closer cluster
//            assignCluster();
//
//            //Calculate new centroids.
//            calculateCentroids();
//
//            iteration++;
//
//            List currentCentroids = getCentroids();
//
//            //Calculates total distance between new and old Centroids
//            double distance = 0;
//            for(int i = 0; i &lt; lastCentroids.size(); i++) {
//                distance += Point.distance(lastCentroids.get(i),currentCentroids.get(i));
//            }
//            System.out.println("#################");
//            System.out.println("Iteration: " + iteration);
//            System.out.println("Centroid distances: " + distance);
//            plotClusters();
//
//            if(distance == 0) {
//                finish = true;
//            }
//        }
//    }
//
//    private void clearClusters() {
//        for(Cluster cluster : clusters) {
//            cluster.clear();
//        }
//    }
//
//    private List getCentroids() {
//        List centroids = new ArrayList(NUM_CLUSTERS);
//        for(Cluster cluster : clusters) {
//            Point aux = cluster.getCentroid();
//            Point point = new Point(aux.getX(),aux.getY());
//            centroids.add(point);
//        }
//        return centroids;
//    }
//
//    private void assignCluster() {
//        double max = Double.MAX_VALUE;
//        double min = max;
//        int cluster = 0;
//        double distance = 0.0;
//
//        for(Point point : points) {
//            min = max;
//            for(int i = 0; i &lt; NUM_CLUSTERS; i++) {
//                Cluster c = clusters.get(i);
//                distance = Point.distance(point, c.getCentroid());
//                if(distance &lt; min){
//                    min = distance;
//                    cluster = i;
//                }
//            }
//            point.setCluster(cluster);
//            clusters.get(cluster).addPoint(point);
//        }
//    }
//
//    private void calculateCentroids() {
//        for(Cluster cluster : clusters) {
//            double sumX = 0;
//            double sumY = 0;
//            List list = cluster.getPoints();
//            int n_points = list.size();
//
//            for(Point point : list) {
//                sumX += point.getX();
//                sumY += point.getY();
//            }
//
//            Point centroid = cluster.getCentroid();
//            if(n_points &gt; 0) {
//                double newX = sumX / n_points;
//                double newY = sumY / n_points;
//                centroid.setX(newX);
//                centroid.setY(newY);
//            }
//        }
//    }
//}