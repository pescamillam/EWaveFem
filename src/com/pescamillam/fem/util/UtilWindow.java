/*
 * PayU Latam - Copyright (c) 2013 - 2017
 * http://www.payulatam.com
 * Date: 21/12/2017
 */
package com.pescamillam.fem.util;

import com.pescamillam.fem.element.Cst;
import com.sun.xml.internal.bind.v2.runtime.reflect.opt.Const;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.util.BigReal;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferStrategy;
import java.math.BigDecimal;
import java.util.List;

/**
 * @author <a href="mailto:p.escamilla@transportsystems.co">Peter Escamilla</a>.
 * @version 0.2.0
 * @since 0.2.0
 */
public class UtilWindow {

	public static void printElements(List<Cst> elements, FieldMatrix<BigReal>[] displacement,
									  FieldMatrix<BigReal>[] speed, FieldMatrix<BigReal>[] acceleration, FieldMatrix<BigReal>[] force) {
		final String title = "Test Window";
		final int width = 30*(Constants.NUM_X+2);
		final int height = 30*(Constants.NUM_Y+3);

		//Creating the frame.
		JFrame frame = new JFrame(title);

		frame.setSize(width, height);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setLocationRelativeTo(null);
		frame.setResizable(false);
		frame.setVisible(true);

		//Creating the canvas.
		Canvas canvas = new Canvas();

		canvas.setSize(width, height);
		canvas.setBackground(Color.BLACK);
		canvas.setVisible(true);
		canvas.setFocusable(false);


		//Putting it all together.
		frame.add(canvas);

		canvas.createBufferStrategy(3);

		boolean running = true;

		BufferStrategy bufferStrategy;
		Graphics graphics;
		int i = 0;

		while (running) {
			bufferStrategy = canvas.getBufferStrategy();
			graphics = bufferStrategy.getDrawGraphics();
			graphics.clearRect(0, 0, width, height);

			graphics.setColor(Color.GREEN);
			graphics.clearRect(0, 0, 30*(Constants.NUM_X+2), 30*(Constants.NUM_Y+3));

			for (int m = 0; m <= Constants.NUM_Y; m++) {
				for (int n = 0; n <= Constants.NUM_X; n++) {
					int x = n*30 + displacement[i].getData()[m*(Constants.NUM_X+1)*2+n*2][0].bigDecimalValue()
							//.multiply(new BigDecimal("500"))
							.intValue();
					int y = m*30 + displacement[i].getData()[m*(Constants.NUM_X+1)*2+n*2+1][0].bigDecimalValue()
							//.multiply(new BigDecimal("500"))
							.intValue();
					graphics.setColor(Color.GREEN);
					graphics.drawOval(x, y, 10, 10);
//                    if (n < numX) {
//                        graphics.drawLine(x + 5, y + 5, x + 35, y + 5);
//                        if (m < numY) {
//                            graphics.drawLine(x + 5, y + 5, x + 35, y + 35);
//                            graphics.drawLine(x + 5, y + 35, x + 35, y + 35);
//                            graphics.drawLine(x + 5, y + 5, x + 5, y + 35);
//                            graphics.drawLine(x + 35, y + 5, x + 35, y + 35);
//                        }
//                    }
					if (speed[i] != null) {
						graphics.setColor(Color.BLUE);
						graphics.drawLine(x + 5, y + 5,
								x + 5 + speed[i].getData()[m*(Constants.NUM_X+1)*2+n*2][0].bigDecimalValue()
										.multiply(new BigDecimal("0.01"))
										.intValue(),
								y + 5 + speed[i].getData()[m*(Constants.NUM_X+1)*2+n*2+1][0].bigDecimalValue()
										.multiply(new BigDecimal("0.01"))
										.intValue());
					}

					if (acceleration[i] != null) {
						graphics.setColor(Color.RED);
						graphics.drawLine(x + 5, y + 5,
								x + 5 + acceleration[i].getData()[m*(Constants.NUM_X+1)*2+n*2][0].bigDecimalValue()
										.multiply(new BigDecimal("0.0001"))
										.intValue(),
								y + 5 + acceleration[i].getData()[m*(Constants.NUM_X+1)*2+n*2+1][0].bigDecimalValue()
										.multiply(new BigDecimal("0.0001"))
										.intValue());
					}

					if (force[i] != null) {
						graphics.setColor(Color.CYAN);
						graphics.drawLine(x + 5, y + 5,
								x + 5 + force[i].getData()[m*(Constants.NUM_X+1)*2+n*2][0].bigDecimalValue()
										.multiply(new BigDecimal("0.001"))
										.intValue(),
								y + 5 + force[i].getData()[m*(Constants.NUM_X+1)*2+n*2+1][0].bigDecimalValue()
										.multiply(new BigDecimal("0.001"))
										.intValue());
					}
				}
			}

			graphics.drawString("t: " + i, 100, 30*(Constants.NUM_Y+1));

			bufferStrategy.show();
			graphics.dispose();
			try {
				Thread.sleep(200L);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			if (i < Constants.NUM_TIMES - 1) {
				i++;
			} else {
				i = 0;
			}
		}
	}
}
